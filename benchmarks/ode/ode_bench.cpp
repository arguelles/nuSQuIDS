// ODE solver benchmark for nuSQuIDS.
//
// Runs one workload with one integration method at one tolerance and reports
// wall time, the number of right-hand-side evaluations, and (optionally) writes
// every flavor expectation value so accuracy can be compared against a
// high-precision reference run.
//
// usage: ode_bench <workload> <method> <tol> [out_file] [reps]
//
// Oscillation parameters: NuFIT 5.2 (2022), normal ordering, without SK
// atmospheric data. Sterile workload: dm2_41 = 1 eV^2, theta_24 = 0.1.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <nuSQuIDS/nuSQuIDS.h>

using namespace nusquids;

static std::atomic<long> g_nfev{0};
static squids::IntegrationStats g_stats;
static long g_nodes = 0;
static double g_evolve_s = 0;  // time spent in EvolveState only (setup excluded)

template <class F>
static void timed(F&& f) {
  auto t0 = std::chrono::steady_clock::now();
  f();
  g_evolve_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}  // energy nodes x neutrino types, for per-node methods

static void add_stats(const squids::IntegrationStats& s) {
  g_stats.rhs_evaluations += s.rhs_evaluations;
  g_stats.accepted_steps += s.accepted_steps;
  g_stats.rejected_steps += s.rejected_steps;
}

// nuSQUIDS calls PreDerive exactly once per right-hand-side evaluation, and
// PreDerive calls AddToPreDerive, so this counts function evaluations.
class CountingNuSQUIDS : public nuSQUIDS {
 public:
  using nuSQUIDS::nuSQUIDS;
  CountingNuSQUIDS() = default;
 protected:
  void AddToPreDerive(double) override { g_nfev.fetch_add(1, std::memory_order_relaxed); }
};

struct Method {
  std::string name;
  const gsl_odeiv2_step_type* step;          // GSL-driven methods
  std::function<void(nuSQUIDS&)> configure;  // extra setup (e.g. non-GSL paths)
  Method() : step(nullptr) {}
  Method(std::string n, const gsl_odeiv2_step_type* s, std::function<void(nuSQUIDS&)> c)
      : name(std::move(n)), step(s), configure(std::move(c)) {}
};

static std::map<std::string, Method> methods() {
  std::map<std::string, Method> m;
  auto add = [&](const std::string& n, const gsl_odeiv2_step_type* s) { m[n] = Method{n, s, nullptr}; };
  add("rkf45", gsl_odeiv2_step_rkf45);
  add("rkck", gsl_odeiv2_step_rkck);
  add("rk8pd", gsl_odeiv2_step_rk8pd);
  add("rk4", gsl_odeiv2_step_rk4);
  add("rk2", gsl_odeiv2_step_rk2);
  add("msadams", gsl_odeiv2_step_msadams);
  // native SQuIDS integrators: "<method>" uses the PI controller,
  // "<method>:I" and "<method>:H211b" select the other controllers
  const std::pair<const char*, squids::ODEMethod> natives[] = {
      {"n_rkf45", squids::ODEMethod::RKF45}, {"n_bs3", squids::ODEMethod::BS3},
      {"n_dp5", squids::ODEMethod::DP5},     {"n_tsit5", squids::ODEMethod::Tsit5},
      {"n_vern6", squids::ODEMethod::Vern6}, {"n_vern7", squids::ODEMethod::Vern7},
      {"n_vern8", squids::ODEMethod::Vern8}, {"n_vern9", squids::ODEMethod::Vern9},
      {"n_rkv76", squids::ODEMethod::RKV76}, {"n_dop853", squids::ODEMethod::DOP853}};
  const std::pair<const char*, squids::StepController> ctls[] = {
      {"", squids::StepController::PI}, {":I", squids::StepController::I},
      {":H211b", squids::StepController::H211b}};
  // exact solution for constant density without interactions (rotation to matter eigenstates)
  m["analytic"] = Method{"analytic", nullptr,
                         [](nuSQUIDS& n) { n.Set_AllowConstantDensityOscillationOnlyEvolution(true); }};
  // per-node exponential integrators (Schrodinger picture)
  m["ni_magnus4"] = Method{"ni_magnus4", nullptr,
                           [](nuSQUIDS& n) { n.Set_NodeEvolution(NodeEvolution::Magnus4); }};
  m["ni_cfmagnus4"] = Method{"ni_cfmagnus4", nullptr,
                             [](nuSQUIDS& n) { n.Set_NodeEvolution(NodeEvolution::CFMagnus4); }};
  m["ni_magnus6"] = Method{"ni_magnus6", nullptr,
                           [](nuSQUIDS& n) { n.Set_NodeEvolution(NodeEvolution::Magnus6); }};
  m["ni_modmagnus"] = Method{"ni_modmagnus", nullptr,
                             [](nuSQUIDS& n) { n.Set_NodeEvolution(NodeEvolution::ModifiedMagnus); }};
  // per-node integration ("ni_" prefix): each energy node and neutrino type separately
  for (const auto& nm : natives)
    for (const auto& c : ctls) {
      const squids::ODEMethod meth = nm.second;
      const squids::StepController ctl = c.second;
      const std::string name = "ni_" + std::string(nm.first + 2) + c.first;
      m[name] = Method{name, nullptr, [meth, ctl](nuSQUIDS& n) {
                         n.Set_ODEMethod(meth);
                         n.Set_StepController(ctl);
                         n.Set_NodeEvolution(NodeEvolution::RungeKutta);
                       }};
    }
  for (const auto& nm : natives)
    for (const auto& c : ctls) {
      const squids::ODEMethod meth = nm.second;
      const squids::StepController ctl = c.second;
      const std::string name = std::string(nm.first) + c.first;
      m[name] = Method{name, nullptr, [meth, ctl](nuSQUIDS& n) {
                         n.Set_ODEMethod(meth);
                         n.Set_StepController(ctl);
                       }};
    }
  return m;
}

static void set_mixing(nuSQUIDS& n, unsigned numneu, const squids::Const& u) {
  n.Set_MixingAngle(0, 1, 33.41 * u.degree);
  n.Set_MixingAngle(0, 2, 8.58 * u.degree);
  n.Set_MixingAngle(1, 2, 42.2 * u.degree);
  n.Set_CPPhase(0, 2, 232. * u.degree);
  n.Set_SquareMassDifference(1, 7.41e-5);
  n.Set_SquareMassDifference(2, 2.507e-3);
  if (numneu > 3) {
    n.Set_SquareMassDifference(3, 1.0);
    n.Set_MixingAngle(1, 3, 0.1);
  }
}

static void configure_solver(nuSQUIDS& n, const Method& m, double tol) {
  n.Set_rel_error(tol);
  n.Set_abs_error(tol);
  if (m.step) n.Set_GSL_step(m.step);
  if (m.configure) m.configure(n);
  n.Set_ProgressBar(false);
}

// Single trajectory through the Earth's diameter (core crossing).
static std::vector<double> run_single(const Method& m, double tol, unsigned numneu,
                                      double Emin_GeV, double Emax_GeV, unsigned ne,
                                      bool interactions, bool constant_density = false) {
  squids::Const u;
  CountingNuSQUIDS n(logspace(Emin_GeV * u.GeV, Emax_GeV * u.GeV, ne), numneu, both, interactions);
  set_mixing(n, numneu, u);
  if (interactions) n.Set_TauRegeneration(true);
  if (constant_density) {
    // control case: no spline knots, integrated numerically (shortcut is off by default)
    n.Set_Body(std::make_shared<ConstantDensity>(5.0, 0.5));
    n.Set_Track(std::make_shared<ConstantDensity::Track>(12742. * u.km));
  } else {
    n.Set_Body(std::make_shared<Earth>());
    n.Set_Track(std::make_shared<Earth::Track>(12742. * u.km));
  }
  configure_solver(n, m, tol);

  marray<double, 3> ini{ne, 2, numneu};
  std::fill(ini.begin(), ini.end(), 0.0);
  for (unsigned ie = 0; ie < ne; ie++)
    for (unsigned r = 0; r < 2; r++) ini[ie][r][1] = 1.0;  // pure nu_mu / nu_mu-bar
  n.Set_initial_state(ini, flavor);
  if (interactions) n.InitializeInteractions();  // cross-section setup is not integration work
  timed([&] { n.EvolveState(); });
  add_stats(n.Get_IntegrationStats());
  g_nodes = ne * 2;

  std::vector<double> out;
  for (unsigned ie = 0; ie < ne; ie++)
    for (unsigned r = 0; r < 2; r++)
      for (unsigned f = 0; f < numneu; f++) out.push_back(n.EvalFlavorAtNode(f, ie, r));
  return out;
}

static std::vector<double> run_atm(const Method& m, double tol, unsigned numneu,
                                   double Emin_GeV, double Emax_GeV, unsigned ne, unsigned ncz,
                                   bool interactions) {
  squids::Const u;
  nuSQUIDSAtm<CountingNuSQUIDS> a(linspace(-1., 0.2, ncz),
                                  logspace(Emin_GeV * u.GeV, Emax_GeV * u.GeV, ne), numneu,
                                  both, interactions);
  a.Set_EvalThreads(1);
  for (unsigned ic = 0; ic < ncz; ic++) {
    auto& n = a.GetnuSQuIDS(ic);
    set_mixing(n, numneu, u);
    if (interactions) n.Set_TauRegeneration(true);
    configure_solver(n, m, tol);
  }
  marray<double, 4> ini{ncz, ne, 2, numneu};
  std::fill(ini.begin(), ini.end(), 0.0);
  for (unsigned ic = 0; ic < ncz; ic++)
    for (unsigned ie = 0; ie < ne; ie++)
      for (unsigned r = 0; r < 2; r++) ini[ic][ie][r][1] = 1.0;
  a.Set_initial_state(ini, flavor);
  timed([&] { a.EvolveState(); });
  for (unsigned ic = 0; ic < ncz; ic++) add_stats(a.GetnuSQuIDS(ic).Get_IntegrationStats());
  g_nodes = ne * 2;

  std::vector<double> out;
  for (unsigned ic = 0; ic < ncz; ic++)
    for (unsigned ie = 0; ie < ne; ie++)
      for (unsigned r = 0; r < 2; r++)
        for (unsigned f = 0; f < numneu; f++)
          out.push_back(a.GetnuSQuIDS(ic).EvalFlavorAtNode(f, ie, r));
  return out;
}

static std::map<std::string, std::function<std::vector<double>(const Method&, double)>> workloads() {
  return {
      // atmospheric-oscillation energies, core-crossing chord
      {"earth_osc", [](const Method& m, double t) { return run_single(m, t, 3, 0.5, 500., 150, false); }},
      // same energies through a constant density: separates spline knots from oscillation
      {"const_osc", [](const Method& m, double t) { return run_single(m, t, 3, 0.5, 500., 150, false, true); }},
      // sub-GeV atmospheric (Hyper-K, ORCA, JUNO): oscillation length below the density scale
      {"earth_subgev", [](const Method& m, double t) { return run_single(m, t, 3, 0.1, 2., 100, false); }},
      // solar day-night (2-20 MeV): oscillation length ~ km
      {"earth_mev", [](const Method& m, double t) { return run_single(m, t, 3, 2.e-3, 2.e-2, 50, false); }},
      // full atmospheric grid, oscillations only
      {"atm_osc", [](const Method& m, double t) { return run_atm(m, t, 3, 1., 1000., 100, 20, false); }},
      // high-energy absorption, NC cascades and tau regeneration through the core
      {"earth_int", [](const Method& m, double t) { return run_single(m, t, 3, 100., 1.e7, 100, true); }},
      // 3+1 sterile, Earth MSW resonance region (IceCube-style)
      {"sterile", [](const Method& m, double t) { return run_single(m, t, 4, 100., 1.e5, 100, false); }},
  };
}

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "usage: ode_bench <workload> <method> <tol> [out_file] [reps]\n";
    return 1;
  }
  const std::string wl = argv[1], mname = argv[2];
  const double tol = std::atof(argv[3]);
  const std::string out_file = argc > 4 ? argv[4] : "";
  const int reps = argc > 5 ? std::atoi(argv[5]) : 1;

  auto ms = methods();
  auto ws = workloads();
  if (!ms.count(mname)) { std::cerr << "unknown method " << mname << "\n"; return 1; }
  if (!ws.count(wl)) { std::cerr << "unknown workload " << wl << "\n"; return 1; }

  std::vector<double> result;
  double best = 1e300;
  long nfev = 0;
  for (int r = 0; r < reps; r++) {
    g_nfev = 0;
    g_stats = squids::IntegrationStats();
    g_evolve_s = 0;
    try {
      result = ws[wl](ms[mname], tol);
    } catch (std::exception& e) {
      std::printf("%s %s %.1e FAILED %s\n", wl.c_str(), mname.c_str(), tol, e.what());
      return 2;
    }
    best = std::min(best, g_evolve_s);
    // per-node methods evaluate one node at a time: report full-system equivalents
    const bool per_node = mname.compare(0, 3, "ni_") == 0;
    nfev = per_node ? static_cast<long>(g_stats.rhs_evaluations / g_nodes) : g_nfev.load();
  }
  std::printf("%s %s %.1e time_s=%.4f nfev=%ld steps=%lu rejected=%lu\n", wl.c_str(), mname.c_str(),
              tol, best, nfev, g_stats.accepted_steps, g_stats.rejected_steps);
  if (!out_file.empty()) {
    FILE* f = std::fopen(out_file.c_str(), "w");
    for (double v : result) std::fprintf(f, "%.17g\n", v);
    std::fclose(f);
  }
  return 0;
}
