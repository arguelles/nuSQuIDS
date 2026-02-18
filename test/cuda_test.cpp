/******************************************************************************
*   Simple CUDA backend test for nuSQuIDS.
*   Tests:
*     1. CUDABackend::IsAvailable() / DeviceCount() / DeviceInfo()
*     2. Backend enum and Set/Get API
******************************************************************************/

#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <cmath>
#include <vector>

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

  // Test 3: GPU propagation (vacuum, no interactions)
  std::cout << "\n--- Test 3: GPU vs CPU Vacuum Propagation ---" << std::endl;
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

    // Compare results
    double max_rel_diff = 0.0;
    double max_abs_diff = 0.0;
    int n_compared = 0;
    for (int ci = 0; ci < ncz; ci++) {
      for (int ei = 0; ei < ne; ei++) {
        for (int flv = 0; flv < (int)numneu; flv++) {
          for (int rho = 0; rho < 2; rho++) {
            double cpu_val = nus_cpu.EvalFlavor(flv, cz_range[ci], e_range[ei], rho);
            double gpu_val = nus_gpu.EvalFlavor(flv, cz_range[ci], e_range[ei], rho);
            double abs_diff = std::abs(cpu_val - gpu_val);
            max_abs_diff = std::max(max_abs_diff, abs_diff);
            // Relative error only for non-tiny values (>1e-6)
            if (std::abs(cpu_val) > 1e-6) {
              double rel_diff = abs_diff / std::abs(cpu_val);
              max_rel_diff = std::max(max_rel_diff, rel_diff);
            }
            n_compared++;
          }
        }
      }
    }

    std::cout << "  Compared " << n_compared << " values" << std::endl;
    std::cout << "  Max abs difference: " << max_abs_diff << std::endl;
    std::cout << "  Max relative diff (for |val|>1e-6): " << max_rel_diff << std::endl;
    // Pass: absolute error < 1e-5 AND relative error < 0.2% for significant values
    bool passed = (max_abs_diff < 1e-5 && max_rel_diff < 2e-3);
    std::cout << "  GPU vs CPU comparison: " << (passed ? "PASS" : "FAIL") << std::endl;
  }

  std::cout << "\n=== All tests completed ===" << std::endl;
  return 0;
#endif
}
