/******************************************************************************
*   CUDA benchmark for nuSQuIDS.
*   Compares GPU vs CPU performance and verifies correctness.
*   Run on a GPU node with: sbatch run_benchmark.sh
******************************************************************************/

#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>

using namespace nusquids;

#ifdef NUSQUIDS_CUDA_ENABLED
#include <nuSQuIDS/cuda/cuda_backend.h>
#endif

struct BenchmarkResult {
  int ncz;
  int ne;
  double cpu_time_ms;
  double gpu_time_ms;
  double speedup;
  double max_rel_diff;
  bool passed;
};

// Helper: setup physics params on a nuSQUIDSAtm
template<typename T>
void setupPhysics(T& nus, unsigned int numneu) {
  nus.Set_MixingAngle(0, 1, 0.563942);
  nus.Set_MixingAngle(0, 2, 0.154085);
  nus.Set_MixingAngle(1, 2, 0.785398);
  nus.Set_SquareMassDifference(1, 7.65e-05);
  nus.Set_SquareMassDifference(2, 0.00247);
  nus.Set_CPPhase(0, 2, 0);
  nus.Set_rel_error(1.0e-6);
  nus.Set_abs_error(1.0e-6);
  nus.Set_IncludeOscillations(true);
}

BenchmarkResult runBenchmark(int ncz, int ne, unsigned int numneu, bool interactions) {
  squids::Const units;
  BenchmarkResult result;
  result.ncz = ncz;
  result.ne = ne;

  double Emin = 1.0e2 * units.GeV;
  double Emax = 1.0e6 * units.GeV;

  auto costh = linspace(-1.0, -0.01, ncz);
  auto energies = logspace(Emin, Emax, ne);

  // Set up initial state: pure mu-neutrino flux
  marray<double,4> inistate{(size_t)ncz, (size_t)ne, 2u, numneu};
  std::fill(inistate.begin(), inistate.end(), 0);
  for (int ci = 0; ci < ncz; ci++)
    for (int ei = 0; ei < ne; ei++)
      for (int rho = 0; rho < 2; rho++)
        inistate[ci][ei][rho][1] = 1.0;

  // ---- CPU evolution ----
  nuSQUIDSAtm<> nus_cpu(costh, energies, numneu, both, interactions);
  setupPhysics(nus_cpu, numneu);
  nus_cpu.Set_initial_state(inistate, flavor);

  auto t0 = std::chrono::high_resolution_clock::now();
  nus_cpu.EvolveState();
  auto t1 = std::chrono::high_resolution_clock::now();
  result.cpu_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

#ifdef NUSQUIDS_CUDA_ENABLED
  // ---- GPU evolution ----
  nuSQUIDSAtm<> nus_gpu(costh, energies, numneu, both, interactions);
  setupPhysics(nus_gpu, numneu);
  nus_gpu.Set_initial_state(inistate, flavor);
  nus_gpu.Set_Backend(Backend::gpu);

  auto t2 = std::chrono::high_resolution_clock::now();
  nus_gpu.EvolveState();
  auto t3 = std::chrono::high_resolution_clock::now();
  result.gpu_time_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

  result.speedup = result.cpu_time_ms / result.gpu_time_ms;

  // ---- Compare results ----
  auto e_range = nus_cpu.GetERange();
  auto cz_range = nus_cpu.GetCosthRange();

  result.max_rel_diff = 0.0;
  double max_abs_diff = 0.0;
  int n_compared = 0;
  int n_failed = 0;
  for (int ci = 0; ci < ncz; ci++) {
    for (int ei = 0; ei < ne; ei++) {
      for (int flv = 0; flv < (int)numneu; flv++) {
        for (int rho = 0; rho < 2; rho++) {
          double cpu_val = nus_cpu.EvalFlavor(flv, cz_range[ci], e_range[ei], rho);
          double gpu_val = nus_gpu.EvalFlavor(flv, cz_range[ci], e_range[ei], rho);
          double abs_diff = std::abs(cpu_val - gpu_val);
          max_abs_diff = std::max(max_abs_diff, abs_diff);
          // Use allclose-style comparison: |diff| / (atol + rtol * max(|a|, |b|))
          // With solver tolerance 1e-6, expect O(1e-4) CPU-GPU differences
          double scale = 1e-4 + 1e-2 * std::max(std::abs(cpu_val), std::abs(gpu_val));
          double scaled_err = abs_diff / scale;
          result.max_rel_diff = std::max(result.max_rel_diff, scaled_err);
          n_compared++;
          if (scaled_err > 1.0) n_failed++;
        }
      }
    }
  }
  // Pass if allclose(atol=1e-4, rtol=1e-2) holds for all values
  result.passed = (n_failed == 0);
#else
  result.gpu_time_ms = 0;
  result.speedup = 0;
  result.max_rel_diff = 0;
  result.passed = false;
#endif

  return result;
}

int main() {
  std::cout << "=== nuSQuIDS CUDA Benchmark ===" << std::endl;

#ifndef NUSQUIDS_CUDA_ENABLED
  std::cout << "CUDA support not compiled in." << std::endl;
  return 0;
#else
  if (!CUDABackend::IsAvailable()) {
    std::cout << "No GPU available." << std::endl;
    return 1;
  }
  std::cout << CUDABackend::DeviceInfo() << std::endl;

  unsigned int numneu = 3;
  bool interactions = false;

  // Benchmark configurations: (ncz, ne)
  struct Config { int ncz; int ne; };
  std::vector<Config> configs = {
    {5, 20},
    {10, 40},
    {20, 60},
    {40, 100},
    {100, 100},
    {100, 200},
    {200, 200},
  };

  std::cout << "\n"
            << std::setw(6) << "ncz"
            << std::setw(6) << "ne"
            << std::setw(14) << "CPU (ms)"
            << std::setw(14) << "GPU (ms)"
            << std::setw(10) << "Speedup"
            << std::setw(14) << "Max Rel Err"
            << std::setw(8) << "Status"
            << std::endl;
  std::cout << std::string(72, '-') << std::endl;

  for (const auto& cfg : configs) {
    std::cout << std::flush;
    auto res = runBenchmark(cfg.ncz, cfg.ne, numneu, interactions);

    std::cout << std::setw(6) << res.ncz
              << std::setw(6) << res.ne
              << std::setw(14) << std::fixed << std::setprecision(1) << res.cpu_time_ms
              << std::setw(14) << std::fixed << std::setprecision(1) << res.gpu_time_ms
              << std::setw(10) << std::fixed << std::setprecision(2) << res.speedup << "x"
              << std::setw(14) << std::scientific << std::setprecision(2) << res.max_rel_diff
              << std::setw(8) << (res.passed ? "PASS" : "FAIL")
              << std::endl;
  }

  std::cout << "\n=== Benchmark complete ===" << std::endl;
  return 0;
#endif
}
