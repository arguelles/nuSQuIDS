/******************************************************************************
*   CUDA interactions test for nuSQuIDS.
*   Compares GPU vs CPU propagation WITH interactions enabled.
*   Tests:
*     1. Atmospheric propagation with CSMS cross-sections (oscillation + NC + CC)
*     2. Same with tau regeneration enabled
*     3. Per-flavor, per-energy diagnostics
******************************************************************************/

#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace nusquids;

#ifdef NUSQUIDS_CUDA_ENABLED
#include <nuSQuIDS/cuda/cuda_backend.h>
#endif

int main() {
  std::cout << "=== nuSQuIDS CUDA Interactions Test ===" << std::endl;

#ifndef NUSQUIDS_CUDA_ENABLED
  std::cout << "CUDA support not compiled in. Skipping." << std::endl;
  return 0;
#else
  if (!CUDABackend::IsAvailable()) {
    std::cout << "No GPU available. Skipping." << std::endl;
    return 0;
  }
  std::cout << "GPU: " << CUDABackend::DeviceInfo() << std::endl;

  squids::Const units;
  const char* flavor_name[] = {"nu_e", "nu_mu", "nu_tau"};
  const char* rho_name[] = {"nu", "nubar"};

  // Test 1: Atmospheric propagation with interactions (no tau regen, no Glashow)
  std::cout << "\n--- Test 1: GPU vs CPU with interactions (NC+CC) ---" << std::endl;
  {
    const unsigned int numneu = 3;
    bool interactions = true;

    double Emin = 1.0e2 * units.GeV;
    double Emax = 1.0e6 * units.GeV;
    double czmin = -1.0;
    double czmax = -0.1;
    int ncz = 5;
    int ne = 40;

    auto costh = linspace(czmin, czmax, ncz);
    auto energies = logspace(Emin, Emax, ne);

    // --- CPU ---
    nuSQUIDSAtm<> nus_cpu(costh, energies, numneu, both, interactions);
    nus_cpu.Set_MixingAngle(0, 1, 0.563942);
    nus_cpu.Set_MixingAngle(0, 2, 0.154085);
    nus_cpu.Set_MixingAngle(1, 2, 0.785398);
    nus_cpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_cpu.Set_SquareMassDifference(2, 0.00247);
    nus_cpu.Set_CPPhase(0, 2, 0);
    nus_cpu.Set_TauRegeneration(false);
    nus_cpu.Set_GlashowResonance(false);

    marray<double,4> inistate{(unsigned int)ncz, (unsigned int)ne, 2u, numneu};
    std::fill(inistate.begin(), inistate.end(), 1);

    nus_cpu.Set_initial_state(inistate, flavor);
    nus_cpu.Set_ProgressBar(false);
    nus_cpu.Set_IncludeOscillations(true);
    nus_cpu.Set_rel_error(1e-6);
    nus_cpu.Set_abs_error(1e-6);

    std::cout << "  Evolving on CPU..." << std::flush;
    nus_cpu.EvolveState();
    std::cout << " done." << std::endl;

    // --- GPU ---
    nuSQUIDSAtm<> nus_gpu(costh, energies, numneu, both, interactions);
    nus_gpu.Set_MixingAngle(0, 1, 0.563942);
    nus_gpu.Set_MixingAngle(0, 2, 0.154085);
    nus_gpu.Set_MixingAngle(1, 2, 0.785398);
    nus_gpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_gpu.Set_SquareMassDifference(2, 0.00247);
    nus_gpu.Set_CPPhase(0, 2, 0);
    nus_gpu.Set_TauRegeneration(false);
    nus_gpu.Set_GlashowResonance(false);

    nus_gpu.Set_initial_state(inistate, flavor);
    nus_gpu.Set_ProgressBar(false);
    nus_gpu.Set_IncludeOscillations(true);
    nus_gpu.Set_rel_error(1e-6);
    nus_gpu.Set_abs_error(1e-6);
    nus_gpu.Set_Backend(Backend::gpu);

    std::cout << "  Evolving on GPU..." << std::flush;
    nus_gpu.EvolveState();
    std::cout << " done." << std::endl;

    // --- Compare ---
    double max_abs_err = 0.0;
    double max_rel_err = 0.0;
    int total_points = 0;
    int oob_count = 0;  // out-of-bounds (flux < 0 or > 1)
    bool pass = true;

    std::cout << std::fixed << std::setprecision(6);

    for (int ic = 0; ic < ncz; ic++) {
      for (int ie = 0; ie < ne; ie++) {
        for (int rho = 0; rho < 2; rho++) {
          for (unsigned int flv = 0; flv < numneu; flv++) {
            double cpu_val = nus_cpu.EvalFlavorAtNode(flv, ic, ie, rho);
            double gpu_val = nus_gpu.EvalFlavorAtNode(flv, ic, ie, rho);

            double abs_err = std::abs(cpu_val - gpu_val);
            double rel_err = (std::abs(cpu_val) > 1e-15) ?
                             abs_err / std::abs(cpu_val) : 0.0;

            max_abs_err = std::max(max_abs_err, abs_err);
            max_rel_err = std::max(max_rel_err, rel_err);
            total_points++;

            if (gpu_val < -1e-6 || gpu_val > 1.0 + 1e-6) oob_count++;

            // Print worst cases
            if (abs_err > 1e-3) {
              std::cout << "  LARGE ERR: cz=" << costh[ic]
                        << " E=" << energies[ie]/units.GeV << "GeV "
                        << rho_name[rho] << " " << flavor_name[flv]
                        << " cpu=" << cpu_val << " gpu=" << gpu_val
                        << " abs=" << abs_err << " rel=" << rel_err
                        << std::endl;
            }
          }
        }
      }
    }

    std::cout << "\n  Summary:" << std::endl;
    std::cout << "    Total comparison points: " << total_points << std::endl;
    std::cout << "    Max absolute error: " << std::scientific << max_abs_err << std::endl;
    std::cout << "    Max relative error: " << max_rel_err << std::endl;
    std::cout << "    Out-of-bounds values: " << oob_count << std::endl;

    // Tolerances for interacting case (relaxed vs oscillation-only)
    double abs_tol = 1e-2;  // interactions introduce larger numerical differences
    double rel_tol = 0.05;  // 5% relative tolerance

    if (max_abs_err > abs_tol || max_rel_err > rel_tol || oob_count > 0) {
      std::cout << "    Result: FAIL (abs_tol=" << abs_tol
                << " rel_tol=" << rel_tol << ")" << std::endl;
      pass = false;
    } else {
      std::cout << "    Result: PASS" << std::endl;
    }

    if (!pass) {
      std::cout << "\n=== INTERACTIONS TEST FAILED ===" << std::endl;
      return 1;
    }
  }

  // Test 2: With tau regeneration
  std::cout << "\n--- Test 2: GPU vs CPU with tau regeneration ---" << std::endl;
  {
    const unsigned int numneu = 3;
    bool interactions = true;

    double Emin = 1.0e4 * units.GeV;
    double Emax = 1.0e8 * units.GeV;
    int ncz = 3;
    int ne = 30;

    auto costh = linspace(-1.0, -0.2, ncz);
    auto energies = logspace(Emin, Emax, ne);

    // --- CPU ---
    nuSQUIDSAtm<> nus_cpu(costh, energies, numneu, both, interactions);
    nus_cpu.Set_MixingAngle(0, 1, 0.563942);
    nus_cpu.Set_MixingAngle(0, 2, 0.154085);
    nus_cpu.Set_MixingAngle(1, 2, 0.785398);
    nus_cpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_cpu.Set_SquareMassDifference(2, 0.00247);
    nus_cpu.Set_CPPhase(0, 2, 0);
    nus_cpu.Set_TauRegeneration(true);
    nus_cpu.Set_GlashowResonance(false);

    marray<double,4> inistate{(unsigned int)ncz, (unsigned int)ne, 2u, numneu};
    std::fill(inistate.begin(), inistate.end(), 0);
    // Inject only tau neutrinos to see tau regen effect
    for (unsigned int ic = 0; ic < (unsigned int)ncz; ic++)
      for (unsigned int ie = 0; ie < (unsigned int)ne; ie++) {
        inistate[ic][ie][0][2] = 1;  // nu_tau
        inistate[ic][ie][1][2] = 1;  // nu_tau_bar
      }

    nus_cpu.Set_initial_state(inistate, flavor);
    nus_cpu.Set_ProgressBar(false);
    nus_cpu.Set_rel_error(1e-6);
    nus_cpu.Set_abs_error(1e-6);

    std::cout << "  Evolving on CPU..." << std::flush;
    nus_cpu.EvolveState();
    std::cout << " done." << std::endl;

    // --- GPU ---
    nuSQUIDSAtm<> nus_gpu(costh, energies, numneu, both, interactions);
    nus_gpu.Set_MixingAngle(0, 1, 0.563942);
    nus_gpu.Set_MixingAngle(0, 2, 0.154085);
    nus_gpu.Set_MixingAngle(1, 2, 0.785398);
    nus_gpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_gpu.Set_SquareMassDifference(2, 0.00247);
    nus_gpu.Set_CPPhase(0, 2, 0);
    nus_gpu.Set_TauRegeneration(true);
    nus_gpu.Set_GlashowResonance(false);

    nus_gpu.Set_initial_state(inistate, flavor);
    nus_gpu.Set_ProgressBar(false);
    nus_gpu.Set_rel_error(1e-6);
    nus_gpu.Set_abs_error(1e-6);
    nus_gpu.Set_Backend(Backend::gpu);

    std::cout << "  Evolving on GPU..." << std::flush;
    nus_gpu.EvolveState();
    std::cout << " done." << std::endl;

    // --- Compare ---
    double max_abs_err = 0.0;
    double max_rel_err = 0.0;
    int total_points = 0;

    for (int ic = 0; ic < ncz; ic++) {
      for (int ie = 0; ie < ne; ie++) {
        for (int rho = 0; rho < 2; rho++) {
          for (unsigned int flv = 0; flv < numneu; flv++) {
            double cpu_val = nus_cpu.EvalFlavorAtNode(flv, ic, ie, rho);
            double gpu_val = nus_gpu.EvalFlavorAtNode(flv, ic, ie, rho);
            double abs_err = std::abs(cpu_val - gpu_val);
            double rel_err = (std::abs(cpu_val) > 1e-15) ?
                             abs_err / std::abs(cpu_val) : 0.0;
            max_abs_err = std::max(max_abs_err, abs_err);
            max_rel_err = std::max(max_rel_err, rel_err);
            total_points++;
          }
        }
      }
    }

    std::cout << "\n  Summary:" << std::endl;
    std::cout << "    Total comparison points: " << total_points << std::endl;
    std::cout << "    Max absolute error: " << std::scientific << max_abs_err << std::endl;
    std::cout << "    Max relative error: " << max_rel_err << std::endl;

    double abs_tol = 5e-2;
    double rel_tol = 0.10;  // 10% for tau regen (more sensitive)
    if (max_abs_err > abs_tol || max_rel_err > rel_tol) {
      std::cout << "    Result: FAIL" << std::endl;
      return 1;
    } else {
      std::cout << "    Result: PASS" << std::endl;
    }
  }

  std::cout << "\n=== ALL INTERACTION TESTS PASSED ===" << std::endl;
  return 0;
#endif
}
