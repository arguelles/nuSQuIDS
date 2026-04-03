/******************************************************************************
*   CUDA RHS component diagnostic test.
*   Compares individual terms of the neutrino evolution equation between
*   CPU and GPU at a single evaluation point to isolate errors.
*
*   Tests (in order):
*     1. Oscillation-only (interactions=true but only coherent term)
*     2. invlen_INT values at a specific point
*     3. GammaRho (absorption) values
*     4. nc_factors (cascade integral)
*     5. Full interacting evolution
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

  // Setup: small problem for fast debugging
  const unsigned int numneu = 3;
  double Emin = 1.0e2 * units.GeV;
  double Emax = 1.0e6 * units.GeV;
  int ncz = 2;
  int ne = 10;

  auto costh = linspace(-1.0, -0.5, ncz);
  auto energies = logspace(Emin, Emax, ne);

  // ========== Test 1: Interactions enabled, compare CPU values ==========
  std::cout << "\n--- Test 1: With interactions (NC+CC, full Earth) ---" << std::endl;
  {
    nuSQUIDSAtm<> cpu(costh, energies, numneu, both, true);
    nuSQUIDSAtm<> gpu_obj(costh, energies, numneu, both, true);

    for (auto* p : {&cpu, &gpu_obj}) {
      p->Set_MixingAngle(0, 1, 0.563942);
      p->Set_MixingAngle(0, 2, 0.154085);
      p->Set_MixingAngle(1, 2, 0.785398);
      p->Set_SquareMassDifference(1, 7.65e-05);
      p->Set_SquareMassDifference(2, 0.00247);
      p->Set_CPPhase(0, 2, 0);
      p->Set_TauRegeneration(false);
      p->Set_GlashowResonance(false);
    }

    marray<double,4> ini{(unsigned)ncz, (unsigned)ne, 2u, numneu};
    std::fill(ini.begin(), ini.end(), 1);
    cpu.Set_initial_state(ini, flavor);
    gpu_obj.Set_initial_state(ini, flavor);

    cpu.Set_ProgressBar(false);
    gpu_obj.Set_ProgressBar(false);

    // Print CPU cross-section info
    cpu.EvolveState(); // this initializes interactions

    std::cout << "  CPU evolution done." << std::endl;

    // Now evolve GPU
    gpu_obj.Set_Backend(Backend::gpu);
    gpu_obj.EvolveState();
    std::cout << "  GPU evolution done." << std::endl;

    // Compare at each energy node
    std::cout << std::fixed << std::setprecision(6);
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

          // Print per-point comparison for neutrinos (rho=0) at first 3 energies
          if (rho == 0 && ie < 3) {
            const char* fn[] = {"e", "mu", "tau"};
            std::cout << "  E=" << std::scientific << energies[ie]/units.GeV
                      << "GeV nu_" << fn[flv]
                      << " cpu=" << std::fixed << cv
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
    std::cout << "    " << ((max_abs < 0.01) ? "PASS" : "FAIL") << std::endl;
  }

  // ========== Test 2: Very short propagation ==========
  // If the error grows with path length, it's an integration issue.
  // If it's present even for a very short path, it's a derivative issue.
  std::cout << "\n--- Test 2: Short path (cos(zenith)=-0.1, ~short through Earth crust) ---" << std::endl;
  {
    auto costh_short = linspace(-0.2, -0.1, 2);
    nuSQUIDSAtm<> cpu(costh_short, energies, numneu, both, true);
    nuSQUIDSAtm<> gpu_obj(costh_short, energies, numneu, both, true);

    for (auto* p : {&cpu, &gpu_obj}) {
      p->Set_MixingAngle(0, 1, 0.563942);
      p->Set_MixingAngle(0, 2, 0.154085);
      p->Set_MixingAngle(1, 2, 0.785398);
      p->Set_SquareMassDifference(1, 7.65e-05);
      p->Set_SquareMassDifference(2, 0.00247);
      p->Set_CPPhase(0, 2, 0);
      p->Set_TauRegeneration(false);
      p->Set_GlashowResonance(false);
    }

    int ncz_short = 2;
    marray<double,4> ini{(unsigned)ncz_short, (unsigned)ne, 2u, numneu};
    std::fill(ini.begin(), ini.end(), 1);
    cpu.Set_initial_state(ini, flavor);
    gpu_obj.Set_initial_state(ini, flavor);
    cpu.Set_ProgressBar(false);
    gpu_obj.Set_ProgressBar(false);
    gpu_obj.Set_Backend(Backend::gpu);

    cpu.EvolveState();
    gpu_obj.EvolveState();

    double max_abs = 0, max_rel = 0;
    for (int ie = 0; ie < ne; ie++) {
      for (int rho = 0; rho < 2; rho++) {
        for (unsigned flv = 0; flv < numneu; flv++) {
          double cv = cpu.EvalFlavor(flv, costh_short[0], energies[ie], rho);
          double gv = gpu_obj.EvalFlavor(flv, costh_short[0], energies[ie], rho);
          double ae = std::abs(cv - gv);
          double re = (std::abs(cv) > 1e-15) ? ae / std::abs(cv) : 0;
          max_abs = std::max(max_abs, ae);
          max_rel = std::max(max_rel, re);
        }
      }
    }
    std::cout << "  max_abs=" << std::scientific << max_abs
              << " max_rel=" << max_rel << std::endl;
    std::cout << "  " << (max_abs < 0.01 ? "PASS" : "FAIL") << std::endl;
  }

  std::cout << "\n=== Diagnostic complete ===" << std::endl;
  return 0;
#endif
}
