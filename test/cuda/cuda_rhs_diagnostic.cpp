/******************************************************************************
*   CUDA RHS component diagnostic test.
*   Compares individual terms of the neutrino evolution equation between
*   CPU and GPU to isolate errors term by term.
*
*   Tests (in order):
*     1. Oscillation-only full evolution — baseline correctness check
*     2. Full interaction evolution — CPU vs GPU nuSQUIDSAtm
*     3. Per-term comparison — DiagnosticNuSQuIDS extracts CPU RHS quantities
******************************************************************************/

#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace nusquids;

#ifdef NUSQUIDS_CUDA_ENABLED
#include <nuSQuIDS/cuda/cuda_backend.h>
#endif

// ===========================================================================
// DiagnosticNuSQuIDS — derived class that exposes protected members for
// term-by-term validation of the RHS.
// ===========================================================================
class DiagnosticNuSQuIDS : public nusquids::nuSQUIDS {
public:
  using nuSQUIDS::nuSQUIDS;

  // Trigger a full PreDerive at position x (sets density, ye, interaction
  // lengths, evolves projectors, etc.)
  void EvaluateAt(double x) { PreDerive(x); }

  // --- Interaction lengths ---
  double get_invlen_INT(int rho, int flv, int ie) const {
    return int_state.invlen_INT[rho][flv][ie];
  }
  double get_invlen_NC(int rho, int flv, int ie) const {
    return int_state.invlen_NC[rho][flv][ie];
  }
  double get_invlen_CC(int rho, int flv, int ie) const {
    return int_state.invlen_CC[rho][flv][ie];
  }

  // --- NC cascade factors ---
  double get_nc_factor(int rho, int flv, int ie) const {
    return nc_factors[rho][flv][ie];
  }

  // --- RHS vector terms ---
  squids::SU_vector get_GammaRho(int ie, int rho) { return GammaRho(ie, rho); }
  squids::SU_vector get_InteractionsRho(int ie, int rho) { return InteractionsRho(ie, rho); }
  squids::SU_vector get_HI(int ie, int rho) { return HI(ie, rho); }

  // --- Evolved projectors ---
  const squids::SU_vector& get_evol_proj(int rho, int flv, int ie) const {
    return evol_b1_proj[rho][flv][ie];
  }

  // --- Body properties at current point ---
  std::vector<double> get_target_ndens() { return GetTargetNumberDensities(); }
  std::vector<double> get_target_fracs() { return GetTargetNumberFractions(); }
  double get_density() const { return current_density; }
  double get_ye() const { return current_ye; }

  // --- Expose for matching ---
  unsigned int get_ne() const { return ne; }
  unsigned int get_numneu() const { return numneu; }
};

// ===========================================================================
// Helpers
// ===========================================================================

void set_osc_params(nuSQUIDSAtm<>& obj) {
  obj.Set_MixingAngle(0, 1, 0.563942);
  obj.Set_MixingAngle(0, 2, 0.154085);
  obj.Set_MixingAngle(1, 2, 0.785398);
  obj.Set_SquareMassDifference(1, 7.65e-05);
  obj.Set_SquareMassDifference(2, 0.00247);
  obj.Set_CPPhase(0, 2, 0);
}

void set_osc_params(DiagnosticNuSQuIDS& obj) {
  obj.Set_MixingAngle(0, 1, 0.563942);
  obj.Set_MixingAngle(0, 2, 0.154085);
  obj.Set_MixingAngle(1, 2, 0.785398);
  obj.Set_SquareMassDifference(1, 7.65e-05);
  obj.Set_SquareMassDifference(2, 0.00247);
  obj.Set_CPPhase(0, 2, 0);
}

// Print SU_vector components (first 9 for SU(3))
void print_su_vector(const std::string& label, const squids::SU_vector& v, int max_comp = 9) {
  auto comps = v.GetComponents();
  int n = std::min((int)comps.size(), max_comp);
  std::cout << "    " << label << ":";
  for (int i = 0; i < n; i++)
    std::cout << " " << std::scientific << std::setprecision(6) << comps[i];
  std::cout << std::endl;
}

int main() {
  std::cout << "=== nuSQuIDS CUDA RHS Component Diagnostic ===" << std::endl;

#ifndef NUSQUIDS_CUDA_ENABLED
  std::cout << "CUDA not compiled in. Skipping." << std::endl;
  return 0;
#else
  if (!CUDABackend::IsAvailable()) {
    std::cout << "No GPU. Skipping." << std::endl;
    return 0;
  }

  squids::Const units;

  // Common parameters — ne=40 avoids alignment issues (rounded_ne==40)
  const unsigned int numneu = 3;
  const double Emin = 1.0e2 * units.GeV;
  const double Emax = 1.0e6 * units.GeV;
  const int ncz = 2;
  const int ne = 40;

  auto costh = linspace(-0.5, -0.1, ncz);
  auto energies = logspace(Emin, Emax, ne);

  // =======================================================================
  // Test 1: Oscillation-only full evolution (interactions=false)
  // Baseline correctness: if oscillation-only fails, the bug is in H0/HI
  // =======================================================================
  std::cout << "\n=== Test 1: Oscillation-only full evolution ===" << std::endl;
  {
    nuSQUIDSAtm<> cpu(costh, energies, numneu, both, false);
    nuSQUIDSAtm<> gpu_obj(costh, energies, numneu, both, false);

    set_osc_params(cpu);
    set_osc_params(gpu_obj);

    marray<double,4> ini{(unsigned)ncz, (unsigned)ne, 2u, numneu};
    std::fill(ini.begin(), ini.end(), 1);
    cpu.Set_initial_state(ini, flavor);
    gpu_obj.Set_initial_state(ini, flavor);

    cpu.Set_ProgressBar(false);
    cpu.Set_rel_error(1e-6);
    cpu.Set_abs_error(1e-6);
    gpu_obj.Set_ProgressBar(false);
    gpu_obj.Set_rel_error(1e-6);
    gpu_obj.Set_abs_error(1e-6);

    cpu.EvolveState();
    std::cout << "  CPU evolution done." << std::endl;

    gpu_obj.Set_Backend(Backend::gpu);
    gpu_obj.EvolveState();
    std::cout << "  GPU evolution done." << std::endl;

    double max_abs = 0, max_rel = 0;
    int oob = 0;
    for (int ie = 0; ie < ne; ie++) {
      for (int rho = 0; rho < 2; rho++) {
        for (unsigned flv = 0; flv < numneu; flv++) {
          double cv = cpu.EvalFlavor(flv, costh[0], energies[ie], rho);
          double gv = gpu_obj.EvalFlavor(flv, costh[0], energies[ie], rho);
          double ae = std::abs(cv - gv);
          double re = (std::abs(cv) > 1e-15) ? ae / std::abs(cv) : 0;
          max_abs = std::max(max_abs, ae);
          max_rel = std::max(max_rel, re);
          if (gv < -1e-6 || gv > 1.0 + 1e-6) oob++;

          if (rho == 0 && ie < 3) {
            const char* fn[] = {"e", "mu", "tau"};
            std::cout << "  E=" << std::scientific << energies[ie]/units.GeV
                      << "GeV nu_" << fn[flv]
                      << " cpu=" << std::fixed << std::setprecision(8) << cv
                      << " gpu=" << gv
                      << " diff=" << std::scientific << ae << std::endl;
          }
        }
      }
    }

    std::cout << "\n  Summary:" << std::endl;
    std::cout << "    max_abs=" << std::scientific << max_abs
              << " max_rel=" << max_rel
              << " oob=" << oob << std::endl;
    bool pass = (max_abs < 1e-5) && (max_rel < 0.002);
    std::cout << "    " << (pass ? "PASS" : "FAIL") << std::endl;
    if (!pass) {
      std::cout << "    WARNING: oscillation-only test failed, interaction tests may be unreliable" << std::endl;
    }
  }

  // =======================================================================
  // Test 2: Full interaction evolution — CPU vs GPU nuSQUIDSAtm
  // interactions=true, tau regen=false, Glashow=false
  // =======================================================================
  std::cout << "\n=== Test 2: Full interaction evolution (NC+CC) ===" << std::endl;
  {
    nuSQUIDSAtm<> cpu(costh, energies, numneu, both, true);
    nuSQUIDSAtm<> gpu_obj(costh, energies, numneu, both, true);

    set_osc_params(cpu);
    set_osc_params(gpu_obj);

    for (auto* p : {&cpu, &gpu_obj}) {
      p->Set_TauRegeneration(false);
      p->Set_GlashowResonance(false);
    }

    marray<double,4> ini{(unsigned)ncz, (unsigned)ne, 2u, numneu};
    std::fill(ini.begin(), ini.end(), 1);
    cpu.Set_initial_state(ini, flavor);
    gpu_obj.Set_initial_state(ini, flavor);

    cpu.Set_ProgressBar(false);
    cpu.Set_rel_error(1e-6);
    cpu.Set_abs_error(1e-6);
    gpu_obj.Set_ProgressBar(false);
    gpu_obj.Set_rel_error(1e-6);
    gpu_obj.Set_abs_error(1e-6);

    cpu.EvolveState();
    std::cout << "  CPU evolution done." << std::endl;

    gpu_obj.Set_Backend(Backend::gpu);
    gpu_obj.EvolveState();
    std::cout << "  GPU evolution done." << std::endl;

    // Compare at first costh bin
    double max_abs = 0, max_rel = 0;
    for (int ie = 0; ie < ne; ie++) {
      for (int rho = 0; rho < 2; rho++) {
        for (unsigned flv = 0; flv < numneu; flv++) {
          double cv = cpu.EvalFlavor(flv, costh[0], energies[ie], rho);
          double gv = gpu_obj.EvalFlavor(flv, costh[0], energies[ie], rho);
          double ae = std::abs(cv - gv);
          double re = (std::abs(cv) > 1e-15) ? ae / std::abs(cv) : 0;
          max_abs = std::max(max_abs, ae);
          max_rel = std::max(max_rel, re);

          if (rho == 0 && ie < 3) {
            const char* fn[] = {"e", "mu", "tau"};
            std::cout << "  E=" << std::scientific << energies[ie]/units.GeV
                      << "GeV nu_" << fn[flv]
                      << " cpu=" << std::fixed << std::setprecision(8) << cv
                      << " gpu=" << gv
                      << " diff=" << std::scientific << ae << std::endl;
          }
        }
      }
    }

    std::cout << "\n  Summary:" << std::endl;
    std::cout << "    max_abs=" << std::scientific << max_abs
              << " max_rel=" << max_rel << std::endl;
    // Fluxes > 1.0 are expected with interactions (NC cascade amplification)
    bool pass = (max_abs < 5e-3) && (max_rel < 0.01);
    std::cout << "    " << (pass ? "PASS" : "FAIL") << std::endl;
  }

  // =======================================================================
  // Test 3: Per-term comparison using DiagnosticNuSQuIDS
  // Sets up a single-baseline nuSQUIDS on the same body/track/initial-state
  // as path 0 of the atmospheric problem, calls PreDerive at midpoint,
  // and prints individual RHS quantities for comparison with GPU.
  // =======================================================================
  std::cout << "\n=== Test 3: Per-term RHS diagnostic (CPU reference values) ===" << std::endl;
  {
    // First, create an atmospheric object to get the EarthAtm body and track
    // for the first zenith bin (costh[0]).
    auto earth = std::make_shared<EarthAtm>();
    auto trk = std::make_shared<EarthAtm::Track>(earth->MakeTrackWithCosine(costh[0]));

    double xini = trk->GetInitialX();
    double xend = trk->GetFinalX();
    double xmid = 0.5 * (xini + xend);

    std::cout << "  Path 0: costh=" << costh[0]
              << " xini=" << std::scientific << xini
              << " xend=" << xend
              << " xmid=" << xmid << std::endl;

    // Create DiagnosticNuSQuIDS for single baseline, same energy grid
    DiagnosticNuSQuIDS diag(energies, numneu, both, true);
    set_osc_params(diag);
    diag.Set_TauRegeneration(false);
    diag.Set_GlashowResonance(false);

    diag.Set_Body(earth);
    diag.Set_Track(trk);

    // Set uniform initial state (all flavors = 1)
    marray<double,3> ini{(unsigned)ne, 2u, numneu};
    std::fill(ini.begin(), ini.end(), 1);
    diag.Set_initial_state(ini, flavor);

    diag.Set_ProgressBar(false);
    diag.Set_rel_error(1e-6);
    diag.Set_abs_error(1e-6);

    // Initialize interactions (fills cross-section tables)
    diag.InitializeInteractions();

    // Now evaluate at midpoint — this calls PreDerive which sets up
    // density, ye, invlen, nc_factors, evolves projectors, etc.
    diag.EvaluateAt(xmid);

    // --- Print body properties ---
    std::cout << "\n  Body at xmid:" << std::endl;
    std::cout << "    density = " << diag.get_density() << " g/cm^3" << std::endl;
    std::cout << "    ye      = " << diag.get_ye() << std::endl;

    auto ndens = diag.get_target_ndens();
    std::cout << "    target number densities:";
    for (size_t t = 0; t < ndens.size(); t++)
      std::cout << " [" << t << "]=" << std::scientific << ndens[t];
    std::cout << std::endl;

    auto fracs = diag.get_target_fracs();
    std::cout << "    target fractions:";
    for (size_t t = 0; t < fracs.size(); t++)
      std::cout << " [" << t << "]=" << fracs[t];
    std::cout << std::endl;

    // --- Print invlen_INT, invlen_CC, invlen_NC ---
    std::cout << "\n  invlen_INT[rho=0][flv][ie] (first 3 energies):" << std::endl;
    for (unsigned flv = 0; flv < numneu; flv++) {
      const char* fn[] = {"e", "mu", "tau"};
      std::cout << "    flv=" << fn[flv] << ":";
      for (int ie = 0; ie < 3 && ie < ne; ie++)
        std::cout << " [" << ie << "]=" << std::scientific << std::setprecision(8)
                  << diag.get_invlen_INT(0, flv, ie);
      std::cout << std::endl;
    }

    std::cout << "\n  invlen_CC[rho=0][flv][ie] (first 3 energies):" << std::endl;
    for (unsigned flv = 0; flv < numneu; flv++) {
      const char* fn[] = {"e", "mu", "tau"};
      std::cout << "    flv=" << fn[flv] << ":";
      for (int ie = 0; ie < 3 && ie < ne; ie++)
        std::cout << " [" << ie << "]=" << std::scientific << std::setprecision(8)
                  << diag.get_invlen_CC(0, flv, ie);
      std::cout << std::endl;
    }

    std::cout << "\n  invlen_NC[rho=0][flv][ie] (first 3 energies):" << std::endl;
    for (unsigned flv = 0; flv < numneu; flv++) {
      const char* fn[] = {"e", "mu", "tau"};
      std::cout << "    flv=" << fn[flv] << ":";
      for (int ie = 0; ie < 3 && ie < ne; ie++)
        std::cout << " [" << ie << "]=" << std::scientific << std::setprecision(8)
                  << diag.get_invlen_NC(0, flv, ie);
      std::cout << std::endl;
    }

    // --- Print nc_factors ---
    std::cout << "\n  nc_factors[rho=0][flv][ie] (first 3 energies):" << std::endl;
    for (unsigned flv = 0; flv < numneu; flv++) {
      const char* fn[] = {"e", "mu", "tau"};
      std::cout << "    flv=" << fn[flv] << ":";
      for (int ie = 0; ie < 3 && ie < ne; ie++)
        std::cout << " [" << ie << "]=" << std::scientific << std::setprecision(8)
                  << diag.get_nc_factor(0, flv, ie);
      std::cout << std::endl;
    }

    // --- Print GammaRho ---
    std::cout << "\n  GammaRho(ie=0, rho=0) components [0..8]:" << std::endl;
    {
      squids::SU_vector gamma = diag.get_GammaRho(0, 0);
      print_su_vector("GammaRho", gamma, 9);
    }

    // --- Print InteractionsRho ---
    std::cout << "\n  InteractionsRho(ie=0, rho=0) components [0..8]:" << std::endl;
    {
      squids::SU_vector inter = diag.get_InteractionsRho(0, 0);
      print_su_vector("InteractionsRho", inter, 9);
    }

    // --- Print HI ---
    std::cout << "\n  HI(ie=0, rho=0) components [0..8]:" << std::endl;
    {
      squids::SU_vector hi = diag.get_HI(0, 0);
      print_su_vector("HI", hi, 9);
    }

    // --- Print a few more energies for GammaRho and InteractionsRho ---
    for (int ie = 1; ie < 3 && ie < ne; ie++) {
      std::cout << "\n  GammaRho(ie=" << ie << ", rho=0) components [0..8]:" << std::endl;
      squids::SU_vector gamma = diag.get_GammaRho(ie, 0);
      print_su_vector("GammaRho", gamma, 9);

      std::cout << "  InteractionsRho(ie=" << ie << ", rho=0) components [0..8]:" << std::endl;
      squids::SU_vector inter = diag.get_InteractionsRho(ie, 0);
      print_su_vector("InteractionsRho", inter, 9);
    }

    // --- Print evolved projectors at ie=0 ---
    std::cout << "\n  evol_b1_proj[rho=0][flv][ie=0] components [0..8]:" << std::endl;
    for (unsigned flv = 0; flv < numneu; flv++) {
      const char* fn[] = {"e", "mu", "tau"};
      squids::SU_vector proj = diag.get_evol_proj(0, flv, 0);
      std::string label = std::string("proj_") + fn[flv];
      print_su_vector(label, proj, 9);
    }

    std::cout << "\n  [These CPU reference values should be compared with GPU kernel printf output]" << std::endl;
  }

  std::cout << "\n=== Diagnostic complete ===" << std::endl;
  return 0;
#endif
}
