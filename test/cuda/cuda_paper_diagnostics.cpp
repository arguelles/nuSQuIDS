/******************************************************************************
*  Per-bin diagnostics + IceCube-style benchmark for the GPU paper.
*  Produces:
*    paper_diag_tau_<jobid>.csv      Test 2 (tau regen) per-bin residuals
*    paper_diag_icecube_<jobid>.csv  IceCube atm grid per-bin residuals
*    paper_diag_timing_<jobid>.txt   wall-clock CPU vs GPU on the IceCube grid
*  Read SLURM_JOB_ID from the env to suffix output files; if unset, uses
*  the current pid so multiple local invocations don't clobber each other.
******************************************************************************/

#include <nuSQuIDS/nuSQuIDS.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unistd.h>

using namespace nusquids;
using clk = std::chrono::steady_clock;

#ifdef NUSQUIDS_CUDA_ENABLED
#include <nuSQuIDS/cuda/cuda_backend.h>
#endif

namespace {

const char* flv_name[] = {"nu_e", "nu_mu", "nu_tau"};
const char* rho_name[] = {"nu", "nubar"};

std::string suffix() {
  if (const char* j = std::getenv("SLURM_JOB_ID")) return std::string(j);
  return std::to_string(static_cast<long>(getpid()));
}

double now_ms(clk::time_point t0) {
  return std::chrono::duration<double, std::milli>(clk::now() - t0).count();
}

// Configure CPU and GPU instances with identical physics.
template <class N>
void common_setup(N& nsq, bool taureg, bool glashow) {
  nsq.Set_MixingAngle(0, 1, 0.563942);
  nsq.Set_MixingAngle(0, 2, 0.154085);
  nsq.Set_MixingAngle(1, 2, 0.785398);
  nsq.Set_SquareMassDifference(1, 7.65e-05);
  nsq.Set_SquareMassDifference(2, 0.00247);
  nsq.Set_CPPhase(0, 2, 0);
  nsq.Set_TauRegeneration(taureg);
  nsq.Set_GlashowResonance(glashow);
  nsq.Set_ProgressBar(false);
  nsq.Set_IncludeOscillations(true);
  nsq.Set_rel_error(1e-6);
  nsq.Set_abs_error(1e-6);
}

// Dump CPU vs GPU EvalFlavor residuals for every (ic, ie, rho, flv) bin.
void dump_residuals(std::ofstream& out,
                    nuSQUIDSAtm<>& cpu, nuSQUIDSAtm<>& gpu,
                    const marray<double, 1>& costh,
                    const marray<double, 1>& energies,
                    int numneu) {
  squids::Const units;
  out << "cz,E_GeV,rho,flavor,cpu,gpu,abs_err,rel_err\n";
  out << std::scientific << std::setprecision(8);
  for (size_t ic = 0; ic < costh.extent(0); ic++) {
    for (size_t ie = 0; ie < energies.extent(0); ie++) {
      for (int rho = 0; rho < 2; rho++) {
        for (int flv = 0; flv < numneu; flv++) {
          double c = cpu.EvalFlavor(flv, costh[ic], energies[ie], rho);
          double g = gpu.EvalFlavor(flv, costh[ic], energies[ie], rho);
          double abs_err = std::abs(c - g);
          double rel_err = (std::abs(c) > 1e-15) ? abs_err / std::abs(c) : 0.0;
          out << costh[ic] << ',' << energies[ie] / units.GeV << ','
              << rho_name[rho] << ',' << flv_name[flv] << ','
              << c << ',' << g << ',' << abs_err << ',' << rel_err << '\n';
        }
      }
    }
  }
}

}  // namespace

int main() {
  std::cout << "=== nuSQuIDS GPU paper diagnostics ===\n";

#ifndef NUSQUIDS_CUDA_ENABLED
  std::cout << "CUDA support not compiled in. Skipping.\n";
  return 0;
#else
  if (!CUDABackend::IsAvailable()) {
    std::cout << "No GPU available. Skipping.\n";
    return 0;
  }
  std::cout << "GPU: " << CUDABackend::DeviceInfo() << "\n";
  squids::Const units;
  const std::string sfx = suffix();

  // ---------------- Part 1: Test 2 (tau regen) per-bin residual dump ---------
  // Mirrors test/cuda/cuda_interactions_test.cpp Test 2 setup so the residual
  // table is directly comparable to the published validation pass/fail.
  {
    std::cout << "\n--- Diagnostic 1: tau regen per-bin dump ---\n";
    const unsigned int numneu = 3;
    const int ncz = 3, ne = 30;
    auto costh = linspace(-1.0, -0.2, ncz);
    auto energies = logspace(1.0e4 * units.GeV, 1.0e8 * units.GeV, ne);

    nuSQUIDSAtm<> cpu(costh, energies, numneu, both, /*interactions=*/true);
    common_setup(cpu, /*taureg=*/true, /*glashow=*/false);
    nuSQUIDSAtm<> gpu(costh, energies, numneu, both, /*interactions=*/true);
    common_setup(gpu, /*taureg=*/true, /*glashow=*/false);

    marray<double, 4> inistate{(unsigned int)ncz, (unsigned int)ne, 2u, numneu};
    std::fill(inistate.begin(), inistate.end(), 0.0);
    for (int ic = 0; ic < ncz; ic++)
      for (int ie = 0; ie < ne; ie++) {
        inistate[ic][ie][0][2] = 1.0;
        inistate[ic][ie][1][2] = 1.0;
      }
    cpu.Set_initial_state(inistate, flavor);
    gpu.Set_initial_state(inistate, flavor);
    gpu.Set_Backend(Backend::gpu);

    std::cout << "  Evolving CPU..." << std::flush;
    cpu.EvolveState();
    std::cout << " done.\n  Evolving GPU..." << std::flush;
    gpu.EvolveState();
    std::cout << " done.\n";

    std::string fn = "paper_diag_tau_" + sfx + ".csv";
    std::ofstream out(fn);
    dump_residuals(out, cpu, gpu, costh, energies, numneu);
    std::cout << "  Wrote " << fn << "\n";
  }

  // ---------------- Part 2: IceCube-style atmospheric benchmark --------------
  // Realistic atmospheric oscillation grid: 20 zenith bins (cosh in [-1, 0]),
  // 200 log-spaced energy bins from 1 GeV to 1 PeV, full interactions.
  // Times CPU vs GPU EvolveState wall-clock, dumps per-bin residuals.
  {
    std::cout << "\n--- Diagnostic 2: IceCube atm benchmark ---\n";
    const unsigned int numneu = 3;
    const int ncz = 20, ne = 200;
    auto costh = linspace(-1.0, -0.05, ncz);   // -0.05 avoids the singular horizon
    auto energies = logspace(1.0 * units.GeV, 1.0e6 * units.GeV, ne);

    nuSQUIDSAtm<> cpu(costh, energies, numneu, both, /*interactions=*/true);
    common_setup(cpu, /*taureg=*/true, /*glashow=*/true);
    nuSQUIDSAtm<> gpu(costh, energies, numneu, both, /*interactions=*/true);
    common_setup(gpu, /*taureg=*/true, /*glashow=*/true);

    // Honda-like atmospheric power-law: E^-2.7 across all flavors / nu types.
    marray<double, 4> inistate{(unsigned int)ncz, (unsigned int)ne, 2u, numneu};
    const double Epivot = 1.0e3 * units.GeV;
    const double gamma = 2.7;
    for (int ic = 0; ic < ncz; ic++)
      for (int ie = 0; ie < ne; ie++) {
        double w = std::pow(energies[ie] / Epivot, -gamma);
        // Atmospheric is dominantly nu_mu / nubar_mu with ~0.5 nu_e / nubar_e.
        for (int rho = 0; rho < 2; rho++) {
          inistate[ic][ie][rho][0] = 0.5 * w;  // nu_e
          inistate[ic][ie][rho][1] = 1.0 * w;  // nu_mu
          inistate[ic][ie][rho][2] = 0.0;      // nu_tau (atmospheric ~0)
        }
      }
    cpu.Set_initial_state(inistate, flavor);
    gpu.Set_initial_state(inistate, flavor);
    gpu.Set_Backend(Backend::gpu);

    auto t0 = clk::now();
    cpu.EvolveState();
    double cpu_ms = now_ms(t0);

    t0 = clk::now();
    gpu.EvolveState();
    double gpu_ms = now_ms(t0);

    std::string ts = "paper_diag_timing_" + sfx + ".txt";
    {
      std::ofstream out(ts);
      out << "ncz=" << ncz << "\n";
      out << "ne=" << ne << "\n";
      out << "Emin_GeV=" << energies[0] / units.GeV << "\n";
      out << "Emax_GeV=" << energies[ne - 1] / units.GeV << "\n";
      out << "spectrum=E^-" << gamma << " (Honda-like, nu_e=0.5, nu_mu=1.0)\n";
      out << "interactions=NC+CC+tau+Glashow\n";
      out << "cpu_ms=" << cpu_ms << "\n";
      out << "gpu_ms=" << gpu_ms << "\n";
      out << "speedup=" << cpu_ms / gpu_ms << "\n";
    }
    std::cout << "  CPU: " << cpu_ms << " ms\n";
    std::cout << "  GPU: " << gpu_ms << " ms (speedup " << cpu_ms / gpu_ms << "x)\n";
    std::cout << "  Wrote " << ts << "\n";

    std::string fn = "paper_diag_icecube_" + sfx + ".csv";
    std::ofstream out(fn);
    dump_residuals(out, cpu, gpu, costh, energies, numneu);
    std::cout << "  Wrote " << fn << "\n";
  }

  // ---------------- Part 3: Oscillation-only IceCube grid ---------------------
  // Same atmospheric grid as Part 2, but with interactions disabled. Isolates
  // the integrator + coherent-evolution cost from the cascade machinery, so
  // the resulting CPU/GPU speedup attributes purely to the oscillation kernel.
  {
    std::cout << "\n--- Diagnostic 3: Oscillation-only IceCube grid ---\n";
    const unsigned int numneu = 3;
    const int ncz = 20, ne = 200;
    auto costh = linspace(-1.0, -0.05, ncz);
    auto energies = logspace(1.0 * units.GeV, 1.0e6 * units.GeV, ne);

    nuSQUIDSAtm<> cpu(costh, energies, numneu, both, /*interactions=*/false);
    common_setup(cpu, /*taureg=*/false, /*glashow=*/false);
    nuSQUIDSAtm<> gpu(costh, energies, numneu, both, /*interactions=*/false);
    common_setup(gpu, /*taureg=*/false, /*glashow=*/false);

    marray<double, 4> inistate{(unsigned int)ncz, (unsigned int)ne, 2u, numneu};
    const double Epivot = 1.0e3 * units.GeV;
    const double gamma = 2.7;
    for (int ic = 0; ic < ncz; ic++)
      for (int ie = 0; ie < ne; ie++) {
        double w = std::pow(energies[ie] / Epivot, -gamma);
        for (int rho = 0; rho < 2; rho++) {
          inistate[ic][ie][rho][0] = 0.5 * w;
          inistate[ic][ie][rho][1] = 1.0 * w;
          inistate[ic][ie][rho][2] = 0.0;
        }
      }
    cpu.Set_initial_state(inistate, flavor);
    gpu.Set_initial_state(inistate, flavor);
    gpu.Set_Backend(Backend::gpu);

    auto t0 = clk::now();
    cpu.EvolveState();
    double cpu_ms = now_ms(t0);

    t0 = clk::now();
    gpu.EvolveState();
    double gpu_ms = now_ms(t0);

    std::string ts = "paper_diag_timing_osc_" + sfx + ".txt";
    {
      std::ofstream out(ts);
      out << "ncz=" << ncz << "\n";
      out << "ne=" << ne << "\n";
      out << "Emin_GeV=" << energies[0] / units.GeV << "\n";
      out << "Emax_GeV=" << energies[ne - 1] / units.GeV << "\n";
      out << "spectrum=E^-" << gamma << " (Honda-like, nu_e=0.5, nu_mu=1.0)\n";
      out << "interactions=DISABLED (oscillation only)\n";
      out << "cpu_ms=" << cpu_ms << "\n";
      out << "gpu_ms=" << gpu_ms << "\n";
      out << "speedup=" << cpu_ms / gpu_ms << "\n";
    }
    std::cout << "  CPU: " << cpu_ms << " ms\n";
    std::cout << "  GPU: " << gpu_ms << " ms (speedup " << cpu_ms / gpu_ms << "x)\n";
    std::cout << "  Wrote " << ts << "\n";

    std::string fn = "paper_diag_icecube_osc_" + sfx + ".csv";
    std::ofstream out(fn);
    dump_residuals(out, cpu, gpu, costh, energies, numneu);
    std::cout << "  Wrote " << fn << "\n";
  }

  std::cout << "\n=== Diagnostics complete ===\n";
  return 0;
#endif
}
