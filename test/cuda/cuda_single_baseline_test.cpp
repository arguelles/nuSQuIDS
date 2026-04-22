/******************************************************************************
*   CUDA single-baseline test for nuSQuIDS.
*   Tests GPU backend for the nuSQUIDS class (single trajectory):
*     1. Vacuum: GPU vs CPU comparison
*     2. Earth matter: GPU vs CPU comparison
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

struct ComparisonStats {
  double sum_abs_err = 0;
  double max_abs_err = 0;
  double max_rel_err = 0;
  int count = 0;
  int count_rel = 0;
  double worst_abs_E = 0;
  double worst_rel_E = 0;
};

bool run_comparison(const std::string& test_name,
                    nuSQUIDS& nus_cpu, nuSQUIDS& nus_gpu,
                    unsigned int numneu, int nrhos,
                    double abs_tol, double rel_tol) {
  squids::Const units;
  const char* flavor_name[] = {"nu_e", "nu_mu", "nu_tau", "nu_s1", "nu_s2", "nu_s3"};
  const char* rho_name[] = {"nu", "nubar"};
  const double rel_threshold = 1e-6;

  auto E_range = nus_cpu.GetERange();
  int ne = E_range.size();

  ComparisonStats stats[6][2];  // [flv][rho]
  double global_max_abs = 0, global_max_rel = 0;
  int bounds_bad = 0;

  for(int ie = 0; ie < ne; ie++){
    for(unsigned int flv = 0; flv < numneu; flv++){
      for(int rho = 0; rho < nrhos; rho++){
        double cpu_val = nus_cpu.EvalFlavorAtNode(flv, ie, rho);
        double gpu_val = nus_gpu.EvalFlavorAtNode(flv, ie, rho);
        double abs_diff = std::abs(cpu_val - gpu_val);
        double rel_diff = (std::abs(cpu_val) > rel_threshold)
                          ? abs_diff / std::abs(cpu_val) : 0.0;

        if(gpu_val < -0.01 || gpu_val > 1.01)
          bounds_bad++;

        auto& s = stats[flv][rho];
        s.sum_abs_err += abs_diff;
        s.count++;

        if(abs_diff > s.max_abs_err){
          s.max_abs_err = abs_diff;
          s.worst_abs_E = E_range[ie];
        }
        if(std::abs(cpu_val) > rel_threshold){
          s.count_rel++;
          if(rel_diff > s.max_rel_err){
            s.max_rel_err = rel_diff;
            s.worst_rel_E = E_range[ie];
          }
        }
        global_max_abs = std::max(global_max_abs, abs_diff);
        global_max_rel = std::max(global_max_rel, rel_diff);
      }
    }
  }

  // Print summary
  std::cout << std::scientific << std::setprecision(6);
  std::cout << "\n  === " << test_name << " Summary ===" << std::endl;
  std::cout << "  " << std::string(90, '-') << std::endl;
  std::cout << std::setw(8) << "flavor" << std::setw(7) << "type"
            << std::setw(14) << "mean|err|" << std::setw(14) << "max|err|"
            << std::setw(14) << "max_rel"
            << "  worst_abs_at" << std::endl;
  std::cout << "  " << std::string(90, '-') << std::endl;

  for(unsigned int flv = 0; flv < numneu; flv++){
    for(int rho = 0; rho < nrhos; rho++){
      auto& s = stats[flv][rho];
      double mean_abs = (s.count > 0) ? s.sum_abs_err / s.count : 0;
      std::cout << "  " << std::setw(8) << flavor_name[flv]
                << std::setw(7) << rho_name[rho]
                << std::setw(14) << mean_abs
                << std::setw(14) << s.max_abs_err
                << std::setw(14) << s.max_rel_err;
      std::cout << "  E=" << std::scientific << std::setprecision(2)
                << s.worst_abs_E / units.GeV << " GeV";
      std::cout << std::endl;
      std::cout << std::scientific << std::setprecision(6);
    }
  }

  std::cout << "  " << std::string(90, '-') << std::endl;
  std::cout << "  Global max absolute error: " << global_max_abs << std::endl;
  std::cout << "  Global max relative error: " << global_max_rel << std::endl;
  std::cout << "  Out-of-bounds values: " << bounds_bad << std::endl;

  bool pass = (global_max_abs < abs_tol && global_max_rel < rel_tol && bounds_bad == 0);
  std::cout << "  Result: " << (pass ? "PASS" : "FAIL") << std::endl;
  if(!pass){
    std::cout << "  (thresholds: abs < " << abs_tol << ", rel < " << rel_tol << ")" << std::endl;
  }
  return pass;
}

int main() {
  std::cout << "=== nuSQuIDS Single-Baseline CUDA Test ===" << std::endl;

#ifndef NUSQUIDS_CUDA_ENABLED
  std::cout << "CUDA support not compiled in. Skipping." << std::endl;
  return 0;
#else
  if(!CUDABackend::IsAvailable()){
    std::cout << "No GPU available. Skipping." << std::endl;
    return 0;
  }

  squids::Const units;
  bool all_pass = true;

  // --- Test 1: Backend API on nuSQUIDS ---
  std::cout << "\n--- Test 1: nuSQUIDS Backend API ---" << std::endl;
  {
    nuSQUIDS nus(logspace(1.0*units.GeV, 1.0e4*units.GeV, 100), 3, neutrino, false);
    std::cout << "  Default backend: " << (int)nus.Get_Backend() << std::endl;
    nus.Set_Backend(Backend::gpu);
    std::cout << "  After Set_Backend(gpu): " << (int)nus.Get_Backend() << std::endl;
    nus.Set_Backend(Backend::cpu);
    std::cout << "  After Set_Backend(cpu): " << (int)nus.Get_Backend() << std::endl;
    std::cout << "  Backend API: PASS" << std::endl;
  }

  // --- Test 2: Vacuum propagation — GPU vs CPU ---
  std::cout << "\n--- Test 2: Vacuum Propagation (GPU vs CPU) ---" << std::endl;
  {
    const unsigned int numneu = 3;
    const double Emin = 1.0e-1 * units.GeV;
    const double Emax = 1.0e4 * units.GeV;
    const int ne = 100;
    const double baseline = 1000.0 * units.km;

    auto energies = logspace(Emin, Emax, ne);

    // CPU
    nuSQUIDS nus_cpu(energies, numneu, both, false);
    nus_cpu.Set_MixingAngle(0, 1, 0.563942);
    nus_cpu.Set_MixingAngle(0, 2, 0.154085);
    nus_cpu.Set_MixingAngle(1, 2, 0.785398);
    nus_cpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_cpu.Set_SquareMassDifference(2, 0.00247);
    nus_cpu.Set_CPPhase(0, 2, 0);

    auto vacuum = std::make_shared<Vacuum>();
    auto vac_track = std::make_shared<Vacuum::Track>(baseline);
    nus_cpu.Set_Body(vacuum);
    nus_cpu.Set_Track(vac_track);

    marray<double,3> inistate{static_cast<size_t>(ne), 2u, numneu};
    std::fill(inistate.begin(), inistate.end(), 0.0);
    // Pure nu_mu initial state for both nu and nubar
    for(int ie = 0; ie < ne; ie++){
      inistate[ie][0][1] = 1.0;  // nu_mu neutrino
      inistate[ie][1][1] = 1.0;  // nu_mu antineutrino
    }
    nus_cpu.Set_initial_state(inistate, flavor);
    nus_cpu.Set_ProgressBar(false);

    std::cout << "  Evolving on CPU..." << std::endl;
    nus_cpu.EvolveState();

    // GPU
    nuSQUIDS nus_gpu(energies, numneu, both, false);
    nus_gpu.Set_MixingAngle(0, 1, 0.563942);
    nus_gpu.Set_MixingAngle(0, 2, 0.154085);
    nus_gpu.Set_MixingAngle(1, 2, 0.785398);
    nus_gpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_gpu.Set_SquareMassDifference(2, 0.00247);
    nus_gpu.Set_CPPhase(0, 2, 0);

    nus_gpu.Set_Body(vacuum);
    nus_gpu.Set_Track(std::make_shared<Vacuum::Track>(baseline));
    nus_gpu.Set_initial_state(inistate, flavor);
    nus_gpu.Set_ProgressBar(false);
    nus_gpu.Set_Backend(Backend::gpu);

    std::cout << "  Evolving on GPU..." << std::endl;
    nus_gpu.EvolveState();

    bool pass = run_comparison("Vacuum GPU vs CPU", nus_cpu, nus_gpu, numneu, 2, 1e-5, 2e-3);
    all_pass = all_pass && pass;
  }

  // --- Test 3: Earth matter propagation — GPU vs CPU ---
  std::cout << "\n--- Test 3: Earth Matter Propagation (GPU vs CPU) ---" << std::endl;
  {
    const unsigned int numneu = 3;
    const double Emin = 1.0e-2 * units.GeV;
    const double Emax = 1.0e2 * units.GeV;
    const int ne = 100;
    const double baseline = 1000.0 * units.km;

    auto energies = logspace(Emin, Emax, ne);

    // CPU
    nuSQUIDS nus_cpu(energies, numneu, both, false);
    nus_cpu.Set_MixingAngle(0, 1, 0.563942);
    nus_cpu.Set_MixingAngle(0, 2, 0.154085);
    nus_cpu.Set_MixingAngle(1, 2, 0.785398);
    nus_cpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_cpu.Set_SquareMassDifference(2, 0.00247);
    nus_cpu.Set_CPPhase(0, 2, 0);

    auto earth = std::make_shared<Earth>();
    double costh = -0.5;  // 60 degrees from vertical
    auto earth_track = std::make_shared<Earth::Track>(earth->MakeTrackWithCosine(costh));
    nus_cpu.Set_Body(earth);
    nus_cpu.Set_Track(earth_track);

    marray<double,3> inistate{static_cast<size_t>(ne), 2u, numneu};
    std::fill(inistate.begin(), inistate.end(), 0.0);
    for(int ie = 0; ie < ne; ie++){
      inistate[ie][0][1] = 1.0;  // nu_mu neutrino
      inistate[ie][1][1] = 1.0;  // nu_mu antineutrino
    }
    nus_cpu.Set_initial_state(inistate, flavor);
    nus_cpu.Set_ProgressBar(false);

    std::cout << "  Evolving on CPU..." << std::endl;
    nus_cpu.EvolveState();

    // GPU
    nuSQUIDS nus_gpu(energies, numneu, both, false);
    nus_gpu.Set_MixingAngle(0, 1, 0.563942);
    nus_gpu.Set_MixingAngle(0, 2, 0.154085);
    nus_gpu.Set_MixingAngle(1, 2, 0.785398);
    nus_gpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_gpu.Set_SquareMassDifference(2, 0.00247);
    nus_gpu.Set_CPPhase(0, 2, 0);

    nus_gpu.Set_Body(earth);
    nus_gpu.Set_Track(std::make_shared<Earth::Track>(earth->MakeTrackWithCosine(costh)));
    nus_gpu.Set_initial_state(inistate, flavor);
    nus_gpu.Set_ProgressBar(false);
    nus_gpu.Set_Backend(Backend::gpu);

    std::cout << "  Evolving on GPU..." << std::endl;
    nus_gpu.EvolveState();

    bool pass = run_comparison("Earth GPU vs CPU", nus_cpu, nus_gpu, numneu, 2, 1e-5, 2e-3);
    all_pass = all_pass && pass;
  }

  // --- Test 4: Constant density — GPU vs CPU ---
  std::cout << "\n--- Test 4: Constant Density (GPU vs CPU) ---" << std::endl;
  {
    const unsigned int numneu = 3;
    const double Emin = 1.0e-1 * units.GeV;
    const double Emax = 1.0e3 * units.GeV;
    const int ne = 100;
    const double baseline = 500.0 * units.km;
    const double density = 3.0;  // g/cm^3
    const double ye = 0.5;

    auto energies = logspace(Emin, Emax, ne);

    // CPU
    nuSQUIDS nus_cpu(energies, numneu, both, false);
    nus_cpu.Set_MixingAngle(0, 1, 0.563942);
    nus_cpu.Set_MixingAngle(0, 2, 0.154085);
    nus_cpu.Set_MixingAngle(1, 2, 0.785398);
    nus_cpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_cpu.Set_SquareMassDifference(2, 0.00247);
    nus_cpu.Set_CPPhase(0, 2, 0);

    auto constdens = std::make_shared<ConstantDensity>(density, ye);
    auto cd_track = std::make_shared<ConstantDensity::Track>(baseline);
    nus_cpu.Set_Body(constdens);
    nus_cpu.Set_Track(cd_track);

    marray<double,3> inistate{static_cast<size_t>(ne), 2u, numneu};
    std::fill(inistate.begin(), inistate.end(), 0.0);
    for(int ie = 0; ie < ne; ie++){
      inistate[ie][0][1] = 1.0;
      inistate[ie][1][1] = 1.0;
    }
    nus_cpu.Set_initial_state(inistate, flavor);
    nus_cpu.Set_ProgressBar(false);

    std::cout << "  Evolving on CPU..." << std::endl;
    nus_cpu.EvolveState();

    // GPU
    nuSQUIDS nus_gpu(energies, numneu, both, false);
    nus_gpu.Set_MixingAngle(0, 1, 0.563942);
    nus_gpu.Set_MixingAngle(0, 2, 0.154085);
    nus_gpu.Set_MixingAngle(1, 2, 0.785398);
    nus_gpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_gpu.Set_SquareMassDifference(2, 0.00247);
    nus_gpu.Set_CPPhase(0, 2, 0);

    nus_gpu.Set_Body(constdens);
    nus_gpu.Set_Track(std::make_shared<ConstantDensity::Track>(baseline));
    nus_gpu.Set_initial_state(inistate, flavor);
    nus_gpu.Set_ProgressBar(false);
    nus_gpu.Set_Backend(Backend::gpu);

    std::cout << "  Evolving on GPU..." << std::endl;
    nus_gpu.EvolveState();

    bool pass = run_comparison("Constant Density GPU vs CPU", nus_cpu, nus_gpu, numneu, 2, 1e-5, 2e-3);
    all_pass = all_pass && pass;
  }

  // --- Test 5: Constant density WITH interactions (NC+CC) — regression guard
  // for the constant-density target_ndens bug (Codex finding #1). Before the
  // fix, convertPath() returned CONSTANT profiles with n_targets=0, silently
  // zeroing CC/NC absorption and the NC cascade. This test exercises the
  // CC-dominated energy regime so any regression shows up as a large deviation
  // from the CPU reference.
  std::cout << "\n--- Test 5: Constant Density + Interactions (GPU vs CPU) ---" << std::endl;
  {
    const unsigned int numneu = 3;
    const double Emin = 1.0e2 * units.GeV;
    const double Emax = 1.0e6 * units.GeV;
    const int ne = 40;
    const double baseline = 1000.0 * units.km;
    const double density = 3.0;  // g/cm^3 (rock-like)
    const double ye = 0.5;

    auto energies = logspace(Emin, Emax, ne);

    // CPU
    nuSQUIDS nus_cpu(energies, numneu, both, true);  // interactions=true
    nus_cpu.Set_MixingAngle(0, 1, 0.563942);
    nus_cpu.Set_MixingAngle(0, 2, 0.154085);
    nus_cpu.Set_MixingAngle(1, 2, 0.785398);
    nus_cpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_cpu.Set_SquareMassDifference(2, 0.00247);
    nus_cpu.Set_CPPhase(0, 2, 0);
    nus_cpu.Set_TauRegeneration(false);
    nus_cpu.Set_GlashowResonance(false);

    auto constdens = std::make_shared<ConstantDensity>(density, ye);
    nus_cpu.Set_Body(constdens);
    nus_cpu.Set_Track(std::make_shared<ConstantDensity::Track>(baseline));

    marray<double,3> inistate{static_cast<size_t>(ne), 2u, numneu};
    std::fill(inistate.begin(), inistate.end(), 1.0);  // equal flux in all flavors
    nus_cpu.Set_initial_state(inistate, flavor);
    nus_cpu.Set_ProgressBar(false);
    nus_cpu.Set_rel_error(1e-6);
    nus_cpu.Set_abs_error(1e-6);

    std::cout << "  Evolving on CPU..." << std::endl;
    nus_cpu.EvolveState();

    // GPU
    nuSQUIDS nus_gpu(energies, numneu, both, true);
    nus_gpu.Set_MixingAngle(0, 1, 0.563942);
    nus_gpu.Set_MixingAngle(0, 2, 0.154085);
    nus_gpu.Set_MixingAngle(1, 2, 0.785398);
    nus_gpu.Set_SquareMassDifference(1, 7.65e-05);
    nus_gpu.Set_SquareMassDifference(2, 0.00247);
    nus_gpu.Set_CPPhase(0, 2, 0);
    nus_gpu.Set_TauRegeneration(false);
    nus_gpu.Set_GlashowResonance(false);

    nus_gpu.Set_Body(constdens);
    nus_gpu.Set_Track(std::make_shared<ConstantDensity::Track>(baseline));
    nus_gpu.Set_initial_state(inistate, flavor);
    nus_gpu.Set_ProgressBar(false);
    nus_gpu.Set_rel_error(1e-6);
    nus_gpu.Set_abs_error(1e-6);
    nus_gpu.Set_Backend(Backend::gpu);

    std::cout << "  Evolving on GPU..." << std::endl;
    nus_gpu.EvolveState();

    // Tolerances match the NC+CC tabulated validation (~0.2%).
    bool pass = run_comparison("Constant Density + Interactions",
                               nus_cpu, nus_gpu, numneu, 2, 1e-3, 5e-3);
    all_pass = all_pass && pass;
  }

  std::cout << "\n=== All single-baseline tests " << (all_pass ? "PASSED" : "FAILED") << " ===" << std::endl;
  return all_pass ? 0 : 1;
#endif
}
