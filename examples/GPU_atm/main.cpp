/******************************************************************************
*   GPU atmospheric example for nuSQuIDS.
*
*   Solves the atmospheric neutrino oscillation problem twice — once on the
*   CPU and once on the CUDA GPU backend — using the same physics inputs, and
*   compares the resulting per-bin fluxes. Prints the maximum absolute and
*   relative disagreement, plus wall-clock time for each side. This is the
*   simplest self-contained illustration of the Set_Backend(Backend::gpu)
*   API added on the cuda-backend branch.
*
*   Build:
*     ./configure --enable-cuda && make examples
*   Run (needs a working CUDA GPU on the host):
*     ./examples/GPU_atm/gpu_atm
*
*   If nuSQuIDS was built without --enable-cuda the GPU section is skipped
*   with a runtime notice and the example still exercises the CPU path.
******************************************************************************/

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include <nuSQuIDS/nuSQuIDS.h>

using namespace nusquids;

namespace {
double flat_flux(double /*enu*/, double /*cz*/) { return 1.0; }
}

int main() {
  squids::Const units;
  const unsigned int numneu = 3;
  const bool interactions = true;  // NC + CC + tau regen + Glashow

  // Grid — kept small enough to run interactively; the paper's benchmark
  // reaches (n_cz, n_E) = (200, 200) at ~50x speedup on a full A100.
  const int ncz = 20;
  const int ne  = 100;
  const double czmin = -1.0;
  const double czmax = -0.1;
  const double emin  = 1.0e2 * units.GeV;
  const double emax  = 1.0e6 * units.GeV;

  auto costh   = linspace(czmin, czmax, ncz);
  auto energies = logspace(emin, emax, ne);

  auto configure = [&](nuSQUIDSAtm<>& nus) {
    nus.Set_MixingAngle(0, 1, 0.563942);
    nus.Set_MixingAngle(0, 2, 0.154085);
    nus.Set_MixingAngle(1, 2, 0.785398);
    nus.Set_SquareMassDifference(1, 7.65e-05);
    nus.Set_SquareMassDifference(2, 0.00247);
    nus.Set_CPPhase(0, 2, 0);

    nus.Set_rel_error(1.0e-6);
    nus.Set_abs_error(1.0e-6);
    nus.Set_ProgressBar(false);
    nus.Set_IncludeOscillations(true);

    marray<double, 4> inistate{(unsigned int)ncz, (unsigned int)ne, 2u, numneu};
    std::fill(inistate.begin(), inistate.end(), 0.0);
    for (int ic = 0; ic < ncz; ic++)
      for (int ie = 0; ie < ne; ie++)
        for (int rho = 0; rho < 2; rho++)
          inistate[ic][ie][rho][1] = flat_flux(energies[ie], costh[ic]);
    nus.Set_initial_state(inistate, flavor);
  };

  auto time_evolve = [](nuSQUIDSAtm<>& nus) {
    auto t0 = std::chrono::steady_clock::now();
    nus.EvolveState();
    auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(t1 - t0).count();
  };

  std::cout << "=== nuSQuIDS GPU atmospheric example ===" << std::endl;
  std::cout << "Grid: " << ncz << " zenith × " << ne
            << " energy, interactions=" << (interactions ? "on" : "off")
            << ", numneu=" << numneu << std::endl;

  // --- CPU reference ---
  nuSQUIDSAtm<> nus_cpu(costh, energies, numneu, both, interactions);
  configure(nus_cpu);
  std::cout << "\nEvolving on CPU..." << std::flush;
  double t_cpu = time_evolve(nus_cpu);
  std::cout << " done in " << t_cpu << " s" << std::endl;

  // --- GPU (skipped cleanly at run time if no CUDA support) ---
  double t_gpu = 0.0;
  bool gpu_ran = false;
  {
    nuSQUIDSAtm<> nus_gpu(costh, energies, numneu, both, interactions);
    configure(nus_gpu);

    try {
      // Set_Backend throws if nuSQuIDS was not built with --enable-cuda,
      // and EvolveState throws if no GPU is available at run time.
      nus_gpu.Set_Backend(Backend::gpu);
      std::cout << "Evolving on GPU..." << std::flush;
      t_gpu = time_evolve(nus_gpu);
      gpu_ran = true;
      std::cout << " done in " << t_gpu << " s" << std::endl;
    } catch (const std::runtime_error& e) {
      std::cout << "\nGPU backend unavailable: " << e.what() << std::endl;
      std::cout << "Skipping GPU vs CPU comparison." << std::endl;
      return 0;
    }

    // --- Compare ---
    double max_abs = 0.0, max_rel = 0.0;
    for (int ic = 0; ic < ncz; ic++) {
      for (int ie = 0; ie < ne; ie++) {
        for (int rho = 0; rho < 2; rho++) {
          for (unsigned int flv = 0; flv < numneu; flv++) {
            double c = nus_cpu.EvalFlavor(flv, costh[ic], energies[ie], rho);
            double g = nus_gpu.EvalFlavor(flv, costh[ic], energies[ie], rho);
            double abs_e = std::abs(c - g);
            double rel_e = (std::abs(c) > 1.0e-15) ? abs_e / std::abs(c) : 0.0;
            max_abs = std::max(max_abs, abs_e);
            max_rel = std::max(max_rel, rel_e);
          }
        }
      }
    }

    std::cout << "\n--- Comparison ---" << std::endl;
    std::cout << "Max absolute error : " << max_abs << std::endl;
    std::cout << "Max relative error : " << max_rel << std::endl;
    if (gpu_ran)
      std::cout << "GPU speedup        : " << (t_cpu / t_gpu) << "×"
                << std::endl;
  }

  return 0;
}
