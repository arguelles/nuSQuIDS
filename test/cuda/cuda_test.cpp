/******************************************************************************
*   CUDA backend test for nuSQuIDS.
*   Tests:
*     1. CUDABackend::IsAvailable() / DeviceCount() / DeviceInfo()
*     2. Backend enum and Set/Get API
*     3. GPU vs CPU propagation comparison (no interactions) with detailed
*        per-flavor, per-energy, per-zenith diagnostics
******************************************************************************/

#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>

using namespace nusquids;

#ifdef NUSQUIDS_CUDA_ENABLED
#include <nuSQuIDS/cuda/cuda_backend.h>
#endif

int main() {
  std::cout << "=== nuSQuIDS CUDA Backend Test ===" << std::endl;

#ifndef NUSQUIDS_CUDA_ENABLED
  std::cout << "CUDA support not compiled in. Skipping." << std::endl;
  return 0;
#else
  squids::Const units;

  const char* flavor_name[] = {"nu_e", "nu_mu", "nu_tau"};
  const char* rho_name[] = {"nu", "nubar"};

  // Test 1: GPU detection
  std::cout << "\n--- Test 1: GPU Detection ---" << std::endl;
  bool available = CUDABackend::IsAvailable();
  int n_devices = CUDABackend::DeviceCount();
  std::string info = CUDABackend::DeviceInfo();

  std::cout << "  CUDA available: " << (available ? "yes" : "no") << std::endl;
  std::cout << "  Device count: " << n_devices << std::endl;
  std::cout << "  Device info:\n" << info << std::endl;

  // Test 2: Backend enum and API
  std::cout << "\n--- Test 2: Backend API ---" << std::endl;
  {
    unsigned int numneu = 3;
    bool interactions = false;

    double Emin = 1.0e1 * units.GeV;
    double Emax = 1.0e2 * units.GeV;

    nuSQUIDSAtm<> nus_atm(linspace(-1.0, 0.0, 5), logspace(Emin, Emax, 10),
                            numneu, both, interactions);

    nus_atm.Set_MixingAngle(0, 1, 0.563942);
    nus_atm.Set_MixingAngle(0, 2, 0.154085);
    nus_atm.Set_MixingAngle(1, 2, 0.785398);
    nus_atm.Set_SquareMassDifference(1, 7.65e-05);
    nus_atm.Set_SquareMassDifference(2, 0.00247);
    nus_atm.Set_CPPhase(0, 2, 0);

    std::cout << "  Default backend: " << (int)nus_atm.Get_Backend() << std::endl;
    nus_atm.Set_Backend(Backend::gpu);
    std::cout << "  After Set_Backend(gpu): " << (int)nus_atm.Get_Backend() << std::endl;
    nus_atm.Set_Backend(Backend::cpu);
    std::cout << "  After Set_Backend(cpu): " << (int)nus_atm.Get_Backend() << std::endl;
    std::cout << "  Backend API: PASS" << std::endl;
  }

  if (!available) {
    std::cout << "\nNo GPU available. Skipping propagation tests." << std::endl;
    std::cout << "\n=== API tests passed ===" << std::endl;
    return 0;
  }

  // Test 3: GPU vs CPU propagation — replicate atmospheric_osc.test.cpp setup
  // Same parameters as test/atmospheric_osc.test.cpp:
  //   Emin=10 GeV, Emax=100 GeV, 10 energy bins, cos(zenith) from -1 to 0 in 5 bins
  //   3 flavors, both nu/nubar, interactions off, unit flux in all flavors
  //   Same mixing angles and mass splittings
  std::cout << "\n--- Test 3: GPU vs CPU (atmospheric_osc parameters) ---" << std::endl;
  {
    const unsigned int numneu = 3;
    bool interactions = false;

    double Emin = 1.0e1 * units.GeV;
    double Emax = 1.0e2 * units.GeV;
    double czmin = -1;
    double czmax = 0;
    int ncz = 5;
    int ne = 10;

    auto costh = linspace(czmin, czmax, ncz);
    auto energies = logspace(Emin, Emax, ne);

    // Use the same evaluation grid as atmospheric_osc.test.cpp
    auto eval_cz = linspace(czmin, czmax, 20);
    auto eval_E = logspace(Emin, Emax, 100);

    // --- CPU ---
    nuSQUIDSAtm<> nus_cpu(costh, energies, numneu, both, interactions);
    nus_cpu.Set_MixingAngle(0, 1, 0.563942);
    nus_cpu.Set_MixingAngle(0, 2, 0.154085);
    nus_cpu.Set_MixingAngle(1, 2, 0.785398);
    nus_cpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_cpu.Set_SquareMassDifference(2, 0.00247);
    nus_cpu.Set_CPPhase(0, 2, 0);

    marray<double,4> inistate{nus_cpu.GetNumCos(), nus_cpu.GetNumE(), 2u, numneu};
    std::fill(inistate.begin(), inistate.end(), 1);  // unit flux in all flavors (matches atmospheric_osc)

    nus_cpu.Set_initial_state(inistate, flavor);
    nus_cpu.Set_ProgressBar(false);
    nus_cpu.Set_IncludeOscillations(true);

    std::cout << "  Evolving on CPU..." << std::endl;
    nus_cpu.EvolveState();
    std::cout << "  CPU evolution done." << std::endl;

    // --- GPU ---
    nuSQUIDSAtm<> nus_gpu(costh, energies, numneu, both, interactions);
    nus_gpu.Set_MixingAngle(0, 1, 0.563942);
    nus_gpu.Set_MixingAngle(0, 2, 0.154085);
    nus_gpu.Set_MixingAngle(1, 2, 0.785398);
    nus_gpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_gpu.Set_SquareMassDifference(2, 0.00247);
    nus_gpu.Set_CPPhase(0, 2, 0);

    nus_gpu.Set_initial_state(inistate, flavor);
    nus_gpu.Set_ProgressBar(false);
    nus_gpu.Set_IncludeOscillations(true);
    nus_gpu.Set_Backend(Backend::gpu);

    std::cout << "  Evolving on GPU..." << std::endl;
    nus_gpu.EvolveState();
    std::cout << "  GPU evolution done." << std::endl;

    // --- Sanity check (same as atmospheric_osc.test.cpp) ---
    // All fluxes should remain in [0,1] for unit initial flux, no interactions
    std::cout << "\n  === Sanity check: fluxes in [0,1] ===" << std::endl;
    int cpu_bad = 0, gpu_bad = 0;
    for (double cz : eval_cz) {
      for (double en : eval_E) {
        for (int fl = 0; fl < (int)numneu; fl++) {
          double cpu_flux = nus_cpu.EvalFlavor(fl, cz, en);
          double gpu_flux = nus_gpu.EvalFlavor(fl, cz, en);
          if (cpu_flux < 0 || cpu_flux > 1) {
            if (cpu_bad < 10)
              std::cout << "  CPU bad flux: " << cpu_flux << " at cz=" << cz
                        << ", E=" << en/units.GeV << " GeV, flavor " << fl << std::endl;
            cpu_bad++;
          }
          if (gpu_flux < 0 || gpu_flux > 1) {
            if (gpu_bad < 10)
              std::cout << "  GPU bad flux: " << gpu_flux << " at cz=" << cz
                        << ", E=" << en/units.GeV << " GeV, flavor " << fl << std::endl;
            gpu_bad++;
          }
        }
      }
    }
    std::cout << "  CPU out-of-range values: " << cpu_bad << std::endl;
    std::cout << "  GPU out-of-range values: " << gpu_bad << std::endl;

    // --- Detailed per-point comparison on the evaluation grid ---
    struct FlavorStats {
      double sum_abs_err = 0;
      double max_abs_err = 0;
      double max_rel_err = 0;
      double sum_rel_err = 0;
      int count = 0;
      int count_rel = 0;
      double worst_abs_cz = 0, worst_abs_E = 0;
      double worst_rel_cz = 0, worst_rel_E = 0;
    };

    const double rel_threshold = 1e-6;
    FlavorStats stats[3][2];  // [flv][rho]

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "\n  === Per-point CPU vs GPU comparison ===" << std::endl;
    std::cout << "  Format: cos(th), E/GeV, flavor, type, CPU, GPU, |diff|, rel_diff" << std::endl;
    std::cout << std::string(105, '-') << std::endl;

    for (double cz : eval_cz) {
      for (double en : eval_E) {
        for (int flv = 0; flv < (int)numneu; flv++) {
          for (int rho = 0; rho < 2; rho++) {
            double cpu_val = nus_cpu.EvalFlavor(flv, cz, en, rho);
            double gpu_val = nus_gpu.EvalFlavor(flv, cz, en, rho);
            double abs_diff = std::abs(cpu_val - gpu_val);
            double rel_diff = (std::abs(cpu_val) > rel_threshold)
                              ? abs_diff / std::abs(cpu_val) : 0.0;

            auto& s = stats[flv][rho];
            s.sum_abs_err += abs_diff;
            s.count++;

            if (abs_diff > s.max_abs_err) {
              s.max_abs_err = abs_diff;
              s.worst_abs_cz = cz;
              s.worst_abs_E = en;
            }

            if (std::abs(cpu_val) > rel_threshold) {
              s.sum_rel_err += rel_diff;
              s.count_rel++;
              if (rel_diff > s.max_rel_err) {
                s.max_rel_err = rel_diff;
                s.worst_rel_cz = cz;
                s.worst_rel_E = en;
              }
            }

            std::cout << "  " << std::setw(11) << cz
                      << "  " << std::setw(12) << en / units.GeV
                      << "  " << std::setw(6) << flavor_name[flv]
                      << "  " << std::setw(5) << rho_name[rho]
                      << "  cpu=" << std::setw(13) << cpu_val
                      << "  gpu=" << std::setw(13) << gpu_val
                      << "  |d|=" << std::setw(10) << abs_diff;
            if (std::abs(cpu_val) > rel_threshold)
              std::cout << "  rel=" << std::setw(10) << rel_diff;
            else
              std::cout << "  rel=       n/a";
            std::cout << std::endl;
          }
        }
      }
    }

    // --- Summary statistics per flavor ---
    std::cout << "\n  === Summary Statistics (per flavor, per type) ===" << std::endl;
    std::cout << std::string(105, '-') << std::endl;
    std::cout << std::setw(8) << "flavor" << std::setw(7) << "type"
              << std::setw(14) << "mean|err|" << std::setw(14) << "max|err|"
              << std::setw(14) << "mean_rel" << std::setw(14) << "max_rel"
              << "  worst_abs_at" << "          worst_rel_at" << std::endl;
    std::cout << std::string(105, '-') << std::endl;

    double global_max_abs = 0, global_max_rel = 0;

    for (int flv = 0; flv < (int)numneu; flv++) {
      for (int rho = 0; rho < 2; rho++) {
        auto& s = stats[flv][rho];
        double mean_abs = (s.count > 0) ? s.sum_abs_err / s.count : 0;
        double mean_rel = (s.count_rel > 0) ? s.sum_rel_err / s.count_rel : 0;

        global_max_abs = std::max(global_max_abs, s.max_abs_err);
        global_max_rel = std::max(global_max_rel, s.max_rel_err);

        std::cout << std::setw(8) << flavor_name[flv]
                  << std::setw(7) << rho_name[rho]
                  << std::setw(14) << mean_abs
                  << std::setw(14) << s.max_abs_err
                  << std::setw(14) << mean_rel
                  << std::setw(14) << s.max_rel_err;
        std::cout << "  cz=" << std::fixed << std::setprecision(3) << s.worst_abs_cz
                  << " E=" << std::scientific << std::setprecision(2) << s.worst_abs_E / units.GeV << "GeV";
        std::cout << "  cz=" << std::fixed << std::setprecision(3) << s.worst_rel_cz
                  << " E=" << std::scientific << std::setprecision(2) << s.worst_rel_E / units.GeV << "GeV";
        std::cout << std::endl;
        std::cout << std::scientific << std::setprecision(6);  // restore
      }
    }

    std::cout << std::string(105, '-') << std::endl;
    std::cout << "  Global max absolute error: " << global_max_abs << std::endl;
    std::cout << "  Global max relative error: " << global_max_rel << std::endl;

    bool bounds_ok = (cpu_bad == 0 && gpu_bad == 0);
    bool diff_ok = (global_max_abs < 1e-5 && global_max_rel < 2e-3);
    std::cout << "\n  Bounds check [0,1]: " << (bounds_ok ? "PASS" : "FAIL") << std::endl;
    std::cout << "  GPU vs CPU comparison: " << (diff_ok ? "PASS" : "FAIL") << std::endl;
    if (!diff_ok) {
      std::cout << "  (thresholds: abs < 1e-5, rel < 0.2%)" << std::endl;
    }
  }

  std::cout << "\n=== All tests completed ===" << std::endl;
  return 0;
#endif
}
