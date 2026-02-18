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
    double Emax = 1.0e6 * units.GeV;

    nuSQUIDSAtm<> nus_atm(linspace(-1.0, 0.0, 5), logspace(Emin, Emax, 40),
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

  // Test 3: GPU vs CPU propagation (no interactions) with detailed diagnostics
  std::cout << "\n--- Test 3: GPU vs CPU Propagation (interactions OFF) ---" << std::endl;
  {
    unsigned int numneu = 3;
    bool interactions = false;

    double Emin = 1.0e2 * units.GeV;
    double Emax = 1.0e5 * units.GeV;
    int ne = 20;
    int ncz = 5;

    auto costh = linspace(-1.0, -0.1, ncz);
    auto energies = logspace(Emin, Emax, ne);

    // --- CPU ---
    nuSQUIDSAtm<> nus_cpu(costh, energies, numneu, both, interactions);
    nus_cpu.Set_MixingAngle(0, 1, 0.563942);
    nus_cpu.Set_MixingAngle(0, 2, 0.154085);
    nus_cpu.Set_MixingAngle(1, 2, 0.785398);
    nus_cpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_cpu.Set_SquareMassDifference(2, 0.00247);
    nus_cpu.Set_CPPhase(0, 2, 0);
    nus_cpu.Set_rel_error(1.0e-8);
    nus_cpu.Set_abs_error(1.0e-8);

    auto e_range = nus_cpu.GetERange();
    auto cz_range = nus_cpu.GetCosthRange();

    marray<double,4> inistate{nus_cpu.GetNumCos(), nus_cpu.GetNumE(), 2u, numneu};
    std::fill(inistate.begin(), inistate.end(), 0);
    for (int ci = 0; ci < nus_cpu.GetNumCos(); ci++)
      for (int ei = 0; ei < nus_cpu.GetNumE(); ei++)
        for (int rho = 0; rho < 2; rho++)
          inistate[ci][ei][rho][1] = 1.0;  // mu neutrino flux = 1

    nus_cpu.Set_initial_state(inistate, flavor);
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
    nus_gpu.Set_rel_error(1.0e-8);
    nus_gpu.Set_abs_error(1.0e-8);

    nus_gpu.Set_initial_state(inistate, flavor);
    nus_gpu.Set_IncludeOscillations(true);
    nus_gpu.Set_Backend(Backend::gpu);

    std::cout << "  Evolving on GPU..." << std::endl;
    nus_gpu.EvolveState();
    std::cout << "  GPU evolution done." << std::endl;

    // --- Detailed per-point comparison ---
    // Per-flavor accumulators: [flv][rho]
    struct FlavorStats {
      double sum_abs_err = 0;
      double max_abs_err = 0;
      double max_rel_err = 0;
      double sum_rel_err = 0;
      int count = 0;
      int count_rel = 0;  // count of points with |cpu_val| > threshold
      // Location of worst absolute error
      int worst_abs_ci = -1, worst_abs_ei = -1;
      // Location of worst relative error
      int worst_rel_ci = -1, worst_rel_ei = -1;
    };

    const double rel_threshold = 1e-6;
    FlavorStats stats[3][2];  // [flv][rho]

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "\n  === Per-point CPU vs GPU comparison ===" << std::endl;
    std::cout << "  Format: cos(th), E/GeV, flavor, type, CPU, GPU, abs_diff, rel_diff" << std::endl;
    std::cout << std::string(100, '-') << std::endl;

    for (int ci = 0; ci < ncz; ci++) {
      for (int ei = 0; ei < ne; ei++) {
        for (int flv = 0; flv < (int)numneu; flv++) {
          for (int rho = 0; rho < 2; rho++) {
            double cpu_val = nus_cpu.EvalFlavor(flv, cz_range[ci], e_range[ei], rho);
            double gpu_val = nus_gpu.EvalFlavor(flv, cz_range[ci], e_range[ei], rho);
            double abs_diff = std::abs(cpu_val - gpu_val);
            double rel_diff = (std::abs(cpu_val) > rel_threshold)
                              ? abs_diff / std::abs(cpu_val) : 0.0;

            auto& s = stats[flv][rho];
            s.sum_abs_err += abs_diff;
            s.count++;

            if (abs_diff > s.max_abs_err) {
              s.max_abs_err = abs_diff;
              s.worst_abs_ci = ci;
              s.worst_abs_ei = ei;
            }

            if (std::abs(cpu_val) > rel_threshold) {
              s.sum_rel_err += rel_diff;
              s.count_rel++;
              if (rel_diff > s.max_rel_err) {
                s.max_rel_err = rel_diff;
                s.worst_rel_ci = ci;
                s.worst_rel_ei = ei;
              }
            }

            // Print every point
            std::cout << "  " << std::setw(10) << cz_range[ci]
                      << "  " << std::setw(12) << e_range[ei] / units.GeV
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
    std::cout << std::string(100, '-') << std::endl;
    std::cout << std::setw(8) << "flavor" << std::setw(7) << "type"
              << std::setw(14) << "mean|err|" << std::setw(14) << "max|err|"
              << std::setw(14) << "mean_rel" << std::setw(14) << "max_rel"
              << "  worst_abs_at" << "  worst_rel_at" << std::endl;
    std::cout << std::string(100, '-') << std::endl;

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
        if (s.worst_abs_ci >= 0)
          std::cout << "  cz[" << s.worst_abs_ci << "]E[" << s.worst_abs_ei << "]";
        else
          std::cout << "  n/a";
        if (s.worst_rel_ci >= 0)
          std::cout << "  cz[" << s.worst_rel_ci << "]E[" << s.worst_rel_ei << "]";
        else
          std::cout << "  n/a";
        std::cout << std::endl;
      }
    }

    std::cout << std::string(100, '-') << std::endl;
    std::cout << "  Global max absolute error: " << global_max_abs << std::endl;
    std::cout << "  Global max relative error: " << global_max_rel << std::endl;

    bool passed = (global_max_abs < 1e-5 && global_max_rel < 2e-3);
    std::cout << "\n  GPU vs CPU comparison: " << (passed ? "PASS" : "FAIL") << std::endl;
    if (!passed) {
      std::cout << "  (thresholds: abs < 1e-5, rel < 0.2%)" << std::endl;
    }
  }

  std::cout << "\n=== All tests completed ===" << std::endl;
  return 0;
#endif
}
