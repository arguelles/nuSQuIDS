// The native integration options of EvolveState reproduce the default (GSL)
// results. Prints nothing when every check passes.
#include <cmath>
#include <iostream>
#include <nuSQuIDS/nuSQuIDS.h>

using namespace nusquids;

static void setup(nuSQUIDS& n, const squids::Const& u){
  n.Set_MixingAngle(0,1,0.5831); n.Set_MixingAngle(0,2,0.1497); n.Set_MixingAngle(1,2,0.7365);
  n.Set_SquareMassDifference(1,7.41e-5); n.Set_SquareMassDifference(2,2.507e-3);
  n.Set_CPPhase(0,2,4.05);
  n.Set_ProgressBar(false);
}

static marray<double,3> muon_state(unsigned int ne){
  marray<double,3> ini{ne,2,3};
  for(unsigned int i = 0; i < ne; i++)
    for(unsigned int r = 0; r < 2; r++)
      for(unsigned int f = 0; f < 3; f++)
        ini[i][r][f] = (f == 1) ? 1.0 : 0.0;
  return ini;
}

// max |P - ref| over energies, neutrino types and flavors
static double run(const marray<double,1>& E, std::shared_ptr<Body> body, std::shared_ptr<Track> track,
                  double tol, void (*configure)(nuSQUIDS&), std::vector<double>* out,
                  const std::vector<double>* ref){
  squids::Const u;
  nuSQUIDS n(E, 3, both, false);
  setup(n, u);
  n.Set_Body(body); n.Set_Track(track);
  n.Set_rel_error(tol); n.Set_abs_error(tol);
  if(configure) configure(n);
  n.Set_initial_state(muon_state(E.extent(0)), flavor);
  n.EvolveState();
  double worst = 0;
  std::vector<double> p;
  for(unsigned int i = 0; i < E.extent(0); i++)
    for(unsigned int r = 0; r < 2; r++)
      for(unsigned int f = 0; f < 3; f++)
        p.push_back(n.EvalFlavorAtNode(f, i, r));
  if(ref)
    for(size_t k = 0; k < p.size(); k++)
      worst = std::max(worst, std::abs(p[k] - (*ref)[k]));
  if(out) *out = p;
  return worst;
}

int main(){
  squids::Const u;
  const marray<double,1> E = logspace(1.*u.GeV, 100.*u.GeV, 30);
  auto earth = std::make_shared<Earth>();
  auto chord = std::make_shared<Earth::Track>(12742.*u.km);

  // Earth reports the PREM features it found: 9 interior radii, each crossed twice
  std::vector<double> bp = earth->GetBreakpoints(*chord);
  if(bp.size() != 18)
    std::cout << "expected 18 breakpoints on a core-crossing chord, got " << bp.size() << std::endl;
  for(size_t i = 0; i < bp.size(); i++)
    if(bp[i] <= 0 || bp[i] >= 12742.*u.km || (i > 0 && bp[i] <= bp[i-1]))
      std::cout << "breakpoint " << i << " out of order or outside the track" << std::endl;

  // reference: the default GSL path at a tight tolerance
  std::vector<double> ref;
  run(E, earth, chord, 1e-12, nullptr, &ref, nullptr);

  struct Case { const char* name; void (*configure)(nuSQUIDS&); double bound; };
  const Case cases[] = {
    {"DP5",  [](nuSQUIDS& n){ n.Set_ODEMethod(squids::ODEMethod::DP5); }, 1e-6},
    {"Tsit5",[](nuSQUIDS& n){ n.Set_ODEMethod(squids::ODEMethod::Tsit5); }, 1e-6},
    {"Vern7",[](nuSQUIDS& n){ n.Set_ODEMethod(squids::ODEMethod::Vern7); }, 1e-6},
    {"Vern9",[](nuSQUIDS& n){ n.Set_ODEMethod(squids::ODEMethod::Vern9); }, 1e-6},
    {"DOP853",[](nuSQUIDS& n){ n.Set_ODEMethod(squids::ODEMethod::DOP853); }, 1e-6},
    {"Vern7 without breakpoints",[](nuSQUIDS& n){ n.Set_ODEMethod(squids::ODEMethod::Vern7); n.Set_UseBreakpoints(false); }, 1e-6},
    {"per-node Vern7",[](nuSQUIDS& n){ n.Set_ODEMethod(squids::ODEMethod::Vern7); n.Set_NodeEvolution(NodeEvolution::RungeKutta); }, 1e-5},
    {"Magnus6",[](nuSQUIDS& n){ n.Set_NodeEvolution(NodeEvolution::Magnus6); }, 1e-5},
    {"modified Magnus",[](nuSQUIDS& n){ n.Set_NodeEvolution(NodeEvolution::ModifiedMagnus); }, 1e-5},
  };
  for(const Case& c : cases){
    double err = run(E, earth, chord, 1e-9, c.configure, nullptr, &ref);
    if(err > c.bound)
      std::cout << c.name << ": max deviation from the GSL reference " << err << std::endl;
  }

  // constant density: the node integrators are exact
  auto cd = std::make_shared<ConstantDensity>(5.0, 0.5);
  auto cdt = std::make_shared<ConstantDensity::Track>(12742.*u.km);
  std::vector<double> exact;
  run(E, cd, cdt, 1e-9, [](nuSQUIDS& n){ n.Set_AllowConstantDensityOscillationOnlyEvolution(true); }, &exact, nullptr);
  const Case exact_cases[] = {
    {"Magnus4",[](nuSQUIDS& n){ n.Set_NodeEvolution(NodeEvolution::Magnus4); }, 1e-12},
    {"Magnus6",[](nuSQUIDS& n){ n.Set_NodeEvolution(NodeEvolution::Magnus6); }, 1e-12},
    {"modified Magnus",[](nuSQUIDS& n){ n.Set_NodeEvolution(NodeEvolution::ModifiedMagnus); }, 1e-12},
  };
  for(const Case& c : exact_cases){
    double err = run(E, cd, cdt, 1e-6, c.configure, nullptr, &exact);
    if(err > c.bound)
      std::cout << c.name << " on constant density: deviation from the exact solution " << err << std::endl;
  }

  // node strategies refuse coupled problems
  {
    nuSQUIDS n(logspace(1.e3*u.GeV, 1.e5*u.GeV, 10), 3, both, true);
    setup(n, u);
    n.Set_Body(earth); n.Set_Track(chord);
    n.Set_NodeEvolution(NodeEvolution::Magnus6);
    n.Set_initial_state(muon_state(10), flavor);
    bool threw = false;
    try{ n.EvolveState(); } catch(std::runtime_error&){ threw = true; }
    if(!threw)
      std::cout << "node evolution with interactions did not throw" << std::endl;
  }

  // time reversal with a native method (breakpoints on a reversed track)
  {
    nuSQUIDS n(logspace(1.e2*u.GeV, 1.e6*u.GeV, 40), 3, both, true);
    setup(n, u);
    auto atm = std::make_shared<EarthAtm>();
    auto tr = std::make_shared<EarthAtm::Track>(atm->MakeTrackWithCosine(-1));
    n.Set_Body(atm); n.Set_Track(tr);
    n.Set_ODEMethod(squids::ODEMethod::Vern7);
    n.Set_rel_error(1e-9); n.Set_abs_error(1e-9);
    marray<double,3> ini = muon_state(40);
    n.Set_initial_state(ini, flavor);
    n.EvolveState();
    tr->ReverseTrack();
    n.EvolveState();
    double worst = 0;
    for(unsigned int i = 0; i < 40; i++)
      for(unsigned int r = 0; r < 2; r++)
        for(unsigned int f = 0; f < 3; f++)
          worst = std::max(worst, std::abs(ini[i][r][f] - n.EvalFlavorAtNode(f, i, r)));
    if(worst > 5e-3)
      std::cout << "time reversal with Vern7: deviation " << worst << std::endl;
  }
  return 0;
}
