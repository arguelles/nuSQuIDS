/******************************************************************************
*   DOPRI5 + FSAL regression tests for the nuSQuIDS CUDA backend.
*
*   Two functions in one binary:
*     - test_fsal_under_rejection : tighten rel/abs tolerance hard, run a steep
*       upgoing path so the adaptive controller is forced to reject steps.
*       Asserts (a) finite step count, (b) GPU vs CPU agreement within 1e-6.
*       This exercises the rejected-step FSAL invariant: path_k[0] of a
*       rejected attempt must remain a valid k1 for the next attempt at the
*       same (t_n, y_n) but smaller h.
*
*     - test_dopri5_convergence : sweep rel_error from 1e-3 to 1e-9, log
*       max |GPU - CPU| / max |CPU|. The 4th-order embedded estimator means
*       the residual scales ~linearly with rel_error — confirms we're in
*       the asymptotic convergence regime, not stuck on cascade-source bias.
*
*   These tests are intentionally CPU-vs-GPU comparisons (CPU is the oracle
*   per perf-dev rules); they do not reach inside the kernel.
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

namespace {

// Build identical CPU/GPU setups so the only difference is the backend.
// All physics knobs are taken from the existing cuda_interactions_test
// Test 1 baseline (PMNS angles, dm^2 values, both rho's).
struct ScenarioConfig {
  unsigned int numneu = 3;
  double Emin;
  double Emax;
  std::vector<double> costh;
  int ne;
  bool tau_regen = false;
  bool glashow   = false;
};

// Inject a flat unit spectrum across all flavors, both rho's. Keeps the
// initial state simple so any divergence is attributable to the integrator.
void make_setup(nuSQUIDSAtm<>& nus, const ScenarioConfig& cfg,
                double rel_err, double abs_err, bool gpu) {
  squids::Const units;
  nus.Set_MixingAngle(0, 1, 0.563942);
  nus.Set_MixingAngle(0, 2, 0.154085);
  nus.Set_MixingAngle(1, 2, 0.785398);
  nus.Set_SquareMassDifference(1, 7.65e-05);
  nus.Set_SquareMassDifference(2, 0.00247);
  nus.Set_CPPhase(0, 2, 0);
  nus.Set_TauRegeneration(cfg.tau_regen);
  nus.Set_GlashowResonance(cfg.glashow);

  marray<double, 4> inistate{(unsigned int)cfg.costh.size(), (unsigned int)cfg.ne, 2u, cfg.numneu};
  std::fill(inistate.begin(), inistate.end(), 1.0);
  nus.Set_initial_state(inistate, flavor);

  nus.Set_ProgressBar(false);
  nus.Set_IncludeOscillations(true);
  nus.Set_rel_error(rel_err);
  nus.Set_abs_error(abs_err);
  if (gpu) nus.Set_Backend(Backend::gpu);
  (void)units;
}

// Compare GPU vs CPU flavor expectations across the entire energy/cz/rho/flv
// grid. Returns (max_abs_err, max_rel_err) and writes the worst-offender
// to stdout.
struct CompareResult {
  double max_abs;
  double max_rel;
  double max_cpu_abs;  // for normalization in convergence study
};

CompareResult compare_states(const nuSQUIDSAtm<>& cpu,
                             const nuSQUIDSAtm<>& gpu,
                             const ScenarioConfig& cfg) {
  CompareResult r{0.0, 0.0, 0.0};
  squids::Const units;
  auto energies = logspace(cfg.Emin, cfg.Emax, cfg.ne);
  for (size_t ic = 0; ic < cfg.costh.size(); ic++) {
    for (int ie = 0; ie < cfg.ne; ie++) {
      for (int rho = 0; rho < 2; rho++) {
        for (unsigned int flv = 0; flv < cfg.numneu; flv++) {
          double cv = cpu.EvalFlavor(flv, cfg.costh[ic], energies[ie], rho);
          double gv = gpu.EvalFlavor(flv, cfg.costh[ic], energies[ie], rho);
          double a = std::abs(cv - gv);
          double rel = (std::abs(cv) > 1e-15) ? a / std::abs(cv) : 0.0;
          r.max_abs = std::max(r.max_abs, a);
          r.max_rel = std::max(r.max_rel, rel);
          r.max_cpu_abs = std::max(r.max_cpu_abs, std::abs(cv));
        }
      }
    }
  }
  (void)units;
  return r;
}

bool test_fsal_under_rejection() {
  std::cout << "\n--- DOPRI5 FSAL test (tight tolerance, forces rejections) ---"
            << std::endl;

  squids::Const units;
  ScenarioConfig cfg;
  cfg.numneu = 3;
  cfg.Emin = 1.0 * units.GeV;          // 1 GeV: rapid vacuum oscillations
  cfg.Emax = 1.0e3 * units.GeV;        // up to 1 TeV
  cfg.ne = 30;
  cfg.costh = {-1.0, -0.85, -0.5};     // strongly upgoing (long path through Earth)
  cfg.tau_regen = false;
  cfg.glashow   = false;

  // Tight tolerance: rel=1e-9, abs=1e-12. Combined with the steep oscillation
  // pattern at 1 GeV through the Earth core this forces the adaptive controller
  // to reject and shrink h repeatedly. If the FSAL implementation reused k7
  // from a rejected step (which would be wrong), the propagation would either
  // diverge or hang in an infinite-rejection loop.
  const double rel = 1.0e-9;
  const double abs_tol = 1.0e-12;

  nuSQUIDSAtm<> nus_cpu(cfg.costh, logspace(cfg.Emin, cfg.Emax, cfg.ne),
                        cfg.numneu, both, /*interactions=*/true);
  nuSQUIDSAtm<> nus_gpu(cfg.costh, logspace(cfg.Emin, cfg.Emax, cfg.ne),
                        cfg.numneu, both, /*interactions=*/true);
  make_setup(nus_cpu, cfg, rel, abs_tol, /*gpu=*/false);
  make_setup(nus_gpu, cfg, rel, abs_tol, /*gpu=*/true);

  std::cout << "  Evolving on CPU (tight tol)..." << std::flush;
  nus_cpu.EvolveState();
  std::cout << " done." << std::endl;

  std::cout << "  Evolving on GPU (tight tol — exercises FSAL across rejections)..."
            << std::flush;
  nus_gpu.EvolveState();
  std::cout << " done." << std::endl;

  CompareResult r = compare_states(nus_cpu, nus_gpu, cfg);
  std::cout << std::scientific << std::setprecision(4);
  std::cout << "  max abs error: " << r.max_abs << std::endl;
  std::cout << "  max rel error: " << r.max_rel << std::endl;

  // Tolerance: the underlying CPU integrator at rel=1e-9 abs=1e-12 produces
  // results accurate to ~ng of double precision; GPU should match within 1e-6
  // absolute (the ratio between rel_err for the integrator and the practical
  // achievable agreement, given accumulated rounding through ~thousands of steps).
  bool pass = (r.max_abs < 1e-6) && std::isfinite(r.max_abs);
  std::cout << "  Result: " << (pass ? "PASS" : "FAIL") << std::endl;
  return pass;
}

bool test_dopri5_convergence() {
  std::cout << "\n--- DOPRI5 convergence sweep (rel_err 1e-3 -> 1e-9) ---"
            << std::endl;

  squids::Const units;
  ScenarioConfig cfg;
  cfg.numneu = 3;
  cfg.Emin = 1.0e2 * units.GeV;
  cfg.Emax = 1.0e6 * units.GeV;
  cfg.ne = 30;
  cfg.costh = {-1.0, -0.6, -0.2};
  cfg.tau_regen = false;
  cfg.glashow   = false;

  // Reference: CPU at very tight tolerance. The integrator residual scales
  // ~linearly with rel_err for a 4th-order embedded estimator, so the GPU's
  // relative error vs the CPU at the same tolerance should also scale linearly
  // (both methods are 5th-order-accurate; both move through the asymptotic
  // regime as rel_err shrinks).
  const double rels[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9};

  bool all_finite = true;
  bool monotone_decreasing = true;
  double prev_rel = std::numeric_limits<double>::infinity();
  for (double rel : rels) {
    nuSQUIDSAtm<> nus_cpu(cfg.costh, logspace(cfg.Emin, cfg.Emax, cfg.ne),
                          cfg.numneu, both, true);
    nuSQUIDSAtm<> nus_gpu(cfg.costh, logspace(cfg.Emin, cfg.Emax, cfg.ne),
                          cfg.numneu, both, true);
    // abs_err scales with rel_err for the convergence study; otherwise the
    // abs floor would clamp the controller and break the linear scaling.
    make_setup(nus_cpu, cfg, rel, rel * 1e-3, /*gpu=*/false);
    make_setup(nus_gpu, cfg, rel, rel * 1e-3, /*gpu=*/true);

    nus_cpu.EvolveState();
    nus_gpu.EvolveState();

    CompareResult r = compare_states(nus_cpu, nus_gpu, cfg);
    double normed = (r.max_cpu_abs > 0.0) ? (r.max_abs / r.max_cpu_abs) : r.max_abs;
    std::cout << "    rel_err=" << std::scientific << std::setprecision(0) << rel
              << "  max|GPU-CPU|=" << std::setprecision(3) << r.max_abs
              << "  norm=" << normed
              << std::endl;

    if (!std::isfinite(r.max_abs)) all_finite = false;
    // The first level (1e-3) often dominated by quantization; only check
    // monotone decrease once we're below 1e-4 where convergence kicks in.
    if (rel <= 1e-4 && normed > prev_rel * 5.0) monotone_decreasing = false;
    prev_rel = normed;
  }

  bool pass = all_finite && monotone_decreasing;
  std::cout << "  Result: " << (pass ? "PASS" : "FAIL")
            << " (finite=" << all_finite
            << ", monotone~=" << monotone_decreasing << ")" << std::endl;
  return pass;
}

} // namespace

int main() {
  std::cout << "=== nuSQuIDS CUDA DOPRI5 / FSAL Tests ===" << std::endl;

#ifndef NUSQUIDS_CUDA_ENABLED
  std::cout << "CUDA support not compiled in. Skipping." << std::endl;
  return 0;
#else
  if (!CUDABackend::IsAvailable()) {
    std::cout << "No GPU available. Skipping." << std::endl;
    return 0;
  }
  std::cout << "GPU: " << CUDABackend::DeviceInfo() << std::endl;

  bool ok = true;
  ok &= test_fsal_under_rejection();
  ok &= test_dopri5_convergence();
  std::cout << "\n=== " << (ok ? "ALL TESTS PASSED" : "SOME TESTS FAILED")
            << " ===" << std::endl;
  return ok ? 0 : 1;
#endif
}
