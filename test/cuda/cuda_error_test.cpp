/******************************************************************************
*   CUDA error-path test for nuSQuIDS.
*   Verifies that host-side preconditions in launchEvolve throw with clear
*   errors instead of silently producing wrong physics.
*     1. SU(4) rejected with std::runtime_error (was silent no-op before 3cf8b80)
*     2. SU(5) rejected with std::runtime_error
*     3. SU(6) rejected with std::runtime_error
*     4. MAX_INT_PAIRS breach rejected with std::runtime_error (bb138c1)
******************************************************************************/

#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <string>
#include <stdexcept>

using namespace nusquids;

#ifdef NUSQUIDS_CUDA_ENABLED
#include <nuSQuIDS/cuda/cuda_backend.h>
#endif

#ifdef NUSQUIDS_CUDA_ENABLED

namespace {

squids::Const g_units;

// Configure a minimal nuSQUIDS instance and try to evolve on the GPU.
// Returns true iff EvolveState() throws std::runtime_error whose message
// contains `expected_fragment`.
bool expect_throw(const std::string& label,
                  unsigned int numneu,
                  int ne,
                  bool interactions,
                  const std::string& expected_fragment) {
  auto energies = logspace(1.0e2 * g_units.GeV, 1.0e6 * g_units.GeV, ne);

  try {
    nuSQUIDS nus(energies, numneu, both, interactions);

    // ConstantDensity gives a non-trivial matter potential and, when
    // interactions=true, populates n_targets>0 so the MAX_INT_PAIRS branch
    // in launchEvolve is actually reached.
    auto body = std::make_shared<ConstantDensity>(5.0, 0.5);
    auto track = std::make_shared<ConstantDensity::Track>(0.0,
                                                          1000.0 * g_units.km);
    nus.Set_Body(body);
    nus.Set_Track(track);

    marray<double, 3> inistate{static_cast<size_t>(ne), 2u, numneu};
    std::fill(inistate.begin(), inistate.end(), 1.0);
    nus.Set_initial_state(inistate, flavor);
    nus.Set_ProgressBar(false);
    nus.Set_TauRegeneration(false);
    nus.Set_GlashowResonance(false);
    nus.Set_rel_error(1e-6);
    nus.Set_abs_error(1e-6);
    nus.Set_Backend(Backend::gpu);

    nus.EvolveState();
  } catch (const std::runtime_error& e) {
    std::string msg = e.what();
    bool ok = msg.find(expected_fragment) != std::string::npos;
    std::cout << "  " << label << ": "
              << (ok ? "PASS" : "FAIL")
              << " (caught: " << msg << ")" << std::endl;
    return ok;
  } catch (const std::exception& e) {
    std::cout << "  " << label
              << ": FAIL (wrong exception type; what()=" << e.what() << ")"
              << std::endl;
    return false;
  }

  std::cout << "  " << label
            << ": FAIL (no exception was thrown; silent no-op regressed)"
            << std::endl;
  return false;
}

}  // namespace

#endif  // NUSQUIDS_CUDA_ENABLED

int main() {
  std::cout << "=== nuSQuIDS CUDA Error-Path Test ===" << std::endl;

#ifndef NUSQUIDS_CUDA_ENABLED
  std::cout << "CUDA support not compiled in. Skipping." << std::endl;
  return 0;
#else
  if (!CUDABackend::IsAvailable()) {
    std::cout << "No GPU available. Skipping." << std::endl;
    return 0;
  }
  std::cout << "GPU: " << CUDABackend::DeviceInfo() << std::endl;

  bool all_pass = true;

  // Tests 1-3: SU(N>3) must fail-fast rather than run a no-op kernel that
  // returns the initial state unchanged. Regression guard for 3cf8b80.
  std::cout << "\n--- SU(N>3) fail-fast (regression guard for 3cf8b80) ---"
            << std::endl;
  all_pass &= expect_throw("Test 1: numneu=4 rejected",
                           /*numneu=*/4, /*ne=*/40,
                           /*interactions=*/false,
                           "SU(4) evolution is not implemented");
  all_pass &= expect_throw("Test 2: numneu=5 rejected",
                           /*numneu=*/5, /*ne=*/40,
                           /*interactions=*/false,
                           "SU(5) evolution is not implemented");
  all_pass &= expect_throw("Test 3: numneu=6 rejected",
                           /*numneu=*/6, /*ne=*/40,
                           /*interactions=*/false,
                           "SU(6) evolution is not implemented");

  // Test 4: MAX_INT_PAIRS host-side check. With EVOLVE_THREADS=128,
  // MAX_INT_PAIRS=4, and nrhos=2 (Type::both), any ne>256 exceeds the bound
  // (2 * ceil(257/128) = 6 > 4). Interactions must be enabled AND
  // n_targets>0 for launchEvolve to reach the check, so ConstantDensity is
  // used above (Vacuum has zero density but still declares targets — either
  // way the check runs). Regression guard for bb138c1.
  std::cout << "\n--- MAX_INT_PAIRS precondition (regression guard for bb138c1) ---"
            << std::endl;
  all_pass &= expect_throw("Test 4: nrhos=2, ne=257 rejected",
                           /*numneu=*/3, /*ne=*/257,
                           /*interactions=*/true,
                           "MAX_INT_PAIRS");

  if (all_pass) {
    std::cout << "\n=== ALL ERROR-PATH TESTS PASSED ===" << std::endl;
    return 0;
  } else {
    std::cout << "\n=== ERROR-PATH TESTS FAILED ===" << std::endl;
    return 1;
  }
#endif
}
