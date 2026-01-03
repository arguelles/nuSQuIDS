/******************************************************************************
*    This program is free software: you can redistribute it and/or modify     *
*   it under the terms of the GNU General Public License as published by      *
*   the Free Software Foundation, either version 3 of the License, or         *
*   (at your option) any later version.                                       *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License         *
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
*                                                                             *
*   Authors:                                                                  *
*      Carlos Arguelles (University of Wisconsin Madison)                     *
*         carguelles@icecube.wisc.edu                                         *
*      Jordi Salvado (University of Wisconsin Madison)                        *
*         jsalvado@icecube.wisc.edu                                           *
*      Christopher Weaver (University of Wisconsin Madison)                   *
*         chris.weaver@icecube.wisc.edu                                       *
******************************************************************************/

/*
 * nuSQuIDS Performance Benchmark Suite
 * =====================================
 *
 * This benchmark tests the performance of nuSQuIDS across different
 * propagation modes and physics configurations:
 *
 * 1. Single Energy Mode
 * 2. Multiple Energy Mode (no interactions)
 * 3. Multiple Energy Mode (with interactions)
 * 4. Atmospheric Mode (no interactions)
 * 5. Atmospheric Mode (with interactions)
 * 6. Atmospheric Mode (with interactions + Glashow)
 * 7. Atmospheric Mode (with interactions + Glashow + Tau Regeneration)
 */

#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>
#include <cmath>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/resources.h>

#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

using namespace nusquids;

// Timing utilities
using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double, std::milli>;

struct BenchmarkResult {
    std::string name;
    std::string description;
    int iterations;
    double total_time_ms;
    double time_per_iter_ms;
    bool passed;
};

// Get system information
std::string get_platform_info() {
#ifdef __APPLE__
    return "macOS";
#elif __linux__
    return "Linux";
#else
    return "Unknown";
#endif
}

std::string get_arch_info() {
#if defined(__aarch64__) || defined(__arm64__)
    return "arm64";
#elif defined(__x86_64__) || defined(_M_X64)
    return "x86_64";
#else
    return "unknown";
#endif
}

// Print utilities
void print_header() {
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "                        nuSQuIDS Performance Benchmark\n";
    std::cout << "                              Version: " << NUSQUIDS_VERSION_STR << "\n";
    std::cout << "================================================================================\n";
    std::cout << "\n";
    std::cout << "System Info:\n";
    std::cout << "  Platform: " << get_platform_info() << " (" << get_arch_info() << ")\n";

    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::cout << "  Date: " << std::ctime(&time);
    std::cout << "\n";
}

void print_test_header(const std::string& name) {
    std::cout << "--------------------------------------------------------------------------------\n";
    std::cout << name << "\n";
    std::cout << "--------------------------------------------------------------------------------\n";
}

void print_test_result(const BenchmarkResult& result) {
    std::cout << "  Test: " << result.description << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "\n";

    if (result.time_per_iter_ms < 1000) {
        std::cout << "  Time per iteration:  " << std::fixed << std::setprecision(2)
                  << result.time_per_iter_ms << " ms\n";
    } else {
        std::cout << "  Time per iteration:  " << std::fixed << std::setprecision(2)
                  << result.time_per_iter_ms / 1000.0 << " s\n";
    }

    if (result.total_time_ms < 1000) {
        std::cout << "  Total time:          " << std::fixed << std::setprecision(0)
                  << result.total_time_ms << " ms\n";
    } else {
        std::cout << "  Total time:          " << std::fixed << std::setprecision(2)
                  << result.total_time_ms / 1000.0 << " s\n";
    }

    std::cout << "  Status: " << (result.passed ? "PASS" : "FAIL") << "\n";
    std::cout << "\n";
}

std::string format_time(double ms) {
    if (ms < 1.0) {
        return std::to_string(static_cast<int>(ms * 1000)) + " us";
    } else if (ms < 1000) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << ms << " ms";
        return oss.str();
    } else {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2) << ms / 1000.0 << " s";
        return oss.str();
    }
}

void print_summary(const std::vector<BenchmarkResult>& results) {
    std::cout << "================================================================================\n";
    std::cout << "                               Summary\n";
    std::cout << "================================================================================\n";

    // Find longest name for alignment
    size_t max_len = 0;
    for (const auto& r : results) {
        max_len = std::max(max_len, r.name.length());
    }

    double total_time = 0;
    int passed = 0;

    for (const auto& r : results) {
        std::cout << "  " << std::left << std::setw(max_len + 2) << r.name;
        std::cout << std::right << std::setw(12) << format_time(r.time_per_iter_ms);
        std::cout << "  " << (r.passed ? "[PASS]" : "[FAIL]") << "\n";
        total_time += r.total_time_ms;
        if (r.passed) passed++;
    }

    std::cout << "\n";
    std::cout << "  Total benchmark time: " << format_time(total_time) << "\n";
    std::cout << "  Tests passed: " << passed << "/" << results.size() << "\n";
    std::cout << "================================================================================\n";
}

// Benchmark functions

BenchmarkResult benchmark_single_energy(int iterations = 100) {
    BenchmarkResult result;
    result.name = "Single Energy (isoscalar)";
    result.description = "Propagate nu_mu through Earth (cos_zen = -1), E = 100 GeV";
    result.iterations = iterations;
    result.passed = true;

    squids::Const units;

    auto start = Clock::now();

    for (int iter = 0; iter < iterations; iter++) {
        nuSQUIDS nus(3, neutrino);

        // Set mixing parameters
        nus.Set_MixingAngle(0, 1, 0.563942);
        nus.Set_MixingAngle(0, 2, 0.154085);
        nus.Set_MixingAngle(1, 2, 0.785398);
        nus.Set_SquareMassDifference(1, 7.65e-05);
        nus.Set_SquareMassDifference(2, 0.00247);
        nus.Set_CPPhase(0, 2, 0.0);

        nus.Set_E(100.0 * units.GeV);

        // Set up Earth atmospheric trajectory
        double phi = acos(-1.0);  // straight up through Earth
        auto earth_atm = std::make_shared<EarthAtm>();
        auto track = std::make_shared<EarthAtm::Track>(earth_atm->MakeTrack(phi));

        nus.Set_Body(earth_atm);
        nus.Set_Track(track);

        // Initial state: pure muon neutrino
        marray<double, 1> ini_state({3}, {0, 1, 0});
        nus.Set_initial_state(ini_state, flavor);

        nus.Set_rel_error(1.0e-10);
        nus.Set_abs_error(1.0e-10);

        nus.EvolveState();

        // Verify result is sensible
        double sum = 0;
        for (int i = 0; i < 3; i++) {
            double val = nus.EvalFlavor(i);
            if (val < -0.01 || val > 1.01) result.passed = false;
            sum += val;
        }
        if (std::abs(sum - 1.0) > 0.01) result.passed = false;
    }

    auto end = Clock::now();
    result.total_time_ms = Duration(end - start).count();
    result.time_per_iter_ms = result.total_time_ms / iterations;

    return result;
}

BenchmarkResult benchmark_multiple_energy_no_interactions(int iterations = 10) {
    BenchmarkResult result;
    result.name = "Multi-E isoscalar (no int.)";
    result.description = "Power-law spectrum through Earth, 1 GeV - 10 TeV (200 nodes)";
    result.iterations = iterations;
    result.passed = true;

    squids::Const units;

    auto start = Clock::now();

    for (int iter = 0; iter < iterations; iter++) {
        nuSQUIDS nus(logspace(1.0 * units.GeV, 1.0e4 * units.GeV, 200), 3, neutrino, false);

        // Set mixing parameters
        nus.Set_MixingAngle(0, 1, 0.563942);
        nus.Set_MixingAngle(0, 2, 0.154085);
        nus.Set_MixingAngle(1, 2, 0.785398);
        nus.Set_SquareMassDifference(1, 7.65e-05);
        nus.Set_SquareMassDifference(2, 0.00247);
        nus.Set_CPPhase(0, 2, 0.0);

        // Set up Earth atmospheric trajectory
        double phi = acos(-1.0);
        auto earth_atm = std::make_shared<EarthAtm>();
        auto track = std::make_shared<EarthAtm::Track>(earth_atm->MakeTrack(phi));

        nus.Set_Body(earth_atm);
        nus.Set_Track(track);

        // Initial state: power-law muon neutrino spectrum
        auto E_range = nus.GetERange();
        marray<double, 2> ini_state{E_range.size(), 3};
        double N0 = 1.0e18;
        for (size_t i = 0; i < E_range.size(); i++) {
            ini_state[i][0] = 0.0;
            ini_state[i][1] = N0 * pow(E_range[i], -2);
            ini_state[i][2] = 0.0;
        }

        nus.Set_initial_state(ini_state, flavor);

        nus.Set_h_max(500.0 * units.km);
        nus.Set_rel_error(1.0e-5);
        nus.Set_abs_error(1.0e-5);
        nus.Set_GSL_step(gsl_odeiv2_step_rk4);

        nus.EvolveState();
    }

    auto end = Clock::now();
    result.total_time_ms = Duration(end - start).count();
    result.time_per_iter_ms = result.total_time_ms / iterations;

    return result;
}

BenchmarkResult benchmark_multiple_energy_with_interactions(int iterations = 5) {
    BenchmarkResult result;
    result.name = "Multi-E isoscalar (with int.)";
    result.description = "Power-law spectrum through Earth, 10 GeV - 1 PeV (100 nodes), NC+CC";
    result.iterations = iterations;
    result.passed = true;

    squids::Const units;

    auto start = Clock::now();

    for (int iter = 0; iter < iterations; iter++) {
        nuSQUIDS nus(logspace(10.0 * units.GeV, 1.0e6 * units.GeV, 100), 3, neutrino, true);

        nus.Set_MixingAngle(0, 1, 0.563942);
        nus.Set_MixingAngle(0, 2, 0.154085);
        nus.Set_MixingAngle(1, 2, 0.785398);
        nus.Set_SquareMassDifference(1, 7.65e-05);
        nus.Set_SquareMassDifference(2, 0.00247);
        nus.Set_CPPhase(0, 2, 0.0);

        double phi = acos(-1.0);
        auto earth_atm = std::make_shared<EarthAtm>();
        auto track = std::make_shared<EarthAtm::Track>(earth_atm->MakeTrack(phi));

        nus.Set_Body(earth_atm);
        nus.Set_Track(track);

        auto E_range = nus.GetERange();
        marray<double, 2> ini_state{E_range.size(), 3};
        double N0 = 1.0e18;
        for (size_t i = 0; i < E_range.size(); i++) {
            ini_state[i][0] = 0.0;
            ini_state[i][1] = N0 * pow(E_range[i], -2);
            ini_state[i][2] = 0.0;
        }

        nus.Set_initial_state(ini_state, flavor);

        nus.Set_h_max(100.0 * units.km);
        nus.Set_rel_error(1.0e-5);
        nus.Set_abs_error(1.0e-5);
        nus.Set_GSL_step(gsl_odeiv2_step_rk4);

        nus.EvolveState();
    }

    auto end = Clock::now();
    result.total_time_ms = Duration(end - start).count();
    result.time_per_iter_ms = result.total_time_ms / iterations;

    return result;
}

BenchmarkResult benchmark_atmospheric_no_interactions(int iterations = 3) {
    BenchmarkResult result;
    result.name = "Atm isoscalar (no int.)";
    result.description = "nuSQUIDSAtm: 10 GeV - 1 PeV (50 E) x 20 zenith, nu+nubar";
    result.iterations = iterations;
    result.passed = true;

    squids::Const units;

    double Emin = 10.0 * units.GeV;
    double Emax = 1.0e6 * units.GeV;
    double czmin = -1.0;
    double czmax = 0.0;

    auto start = Clock::now();

    for (int iter = 0; iter < iterations; iter++) {
        nuSQUIDSAtm<> nus_atm(linspace(czmin, czmax, 20), logspace(Emin, Emax, 50), 3, both, false);

        nus_atm.Set_MixingAngle(0, 1, 0.563942);
        nus_atm.Set_MixingAngle(0, 2, 0.154085);
        nus_atm.Set_MixingAngle(1, 2, 0.785398);
        nus_atm.Set_SquareMassDifference(1, 7.65e-05);
        nus_atm.Set_SquareMassDifference(2, 0.00247);
        nus_atm.Set_CPPhase(0, 2, 0.0);

        nus_atm.Set_rel_error(1.0e-6);
        nus_atm.Set_abs_error(1.0e-6);
        nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);

        // Initial state: flat muon neutrino flux
        marray<double, 4> ini_state{nus_atm.GetNumCos(), nus_atm.GetNumE(), 2, 3};
        std::fill(ini_state.begin(), ini_state.end(), 0.0);
        for (size_t ci = 0; ci < nus_atm.GetNumCos(); ci++) {
            for (size_t ei = 0; ei < nus_atm.GetNumE(); ei++) {
                for (int rho = 0; rho < 2; rho++) {
                    ini_state[ci][ei][rho][1] = 1.0;  // muon flavor
                }
            }
        }

        nus_atm.Set_initial_state(ini_state, flavor);
        nus_atm.Set_IncludeOscillations(true);

        nus_atm.EvolveState();
    }

    auto end = Clock::now();
    result.total_time_ms = Duration(end - start).count();
    result.time_per_iter_ms = result.total_time_ms / iterations;

    return result;
}

BenchmarkResult benchmark_atmospheric_with_interactions(int iterations = 2) {
    BenchmarkResult result;
    result.name = "Atm isoscalar (with int.)";
    result.description = "nuSQUIDSAtm: 10 GeV - 1 PeV (50 E) x 20 zenith, nu+nubar, NC+CC";
    result.iterations = iterations;
    result.passed = true;

    squids::Const units;

    double Emin = 10.0 * units.GeV;
    double Emax = 1.0e6 * units.GeV;
    double czmin = -1.0;
    double czmax = 0.0;

    auto start = Clock::now();

    for (int iter = 0; iter < iterations; iter++) {
        nuSQUIDSAtm<> nus_atm(linspace(czmin, czmax, 20), logspace(Emin, Emax, 50), 3, both, true);

        nus_atm.Set_MixingAngle(0, 1, 0.563942);
        nus_atm.Set_MixingAngle(0, 2, 0.154085);
        nus_atm.Set_MixingAngle(1, 2, 0.785398);
        nus_atm.Set_SquareMassDifference(1, 7.65e-05);
        nus_atm.Set_SquareMassDifference(2, 0.00247);
        nus_atm.Set_CPPhase(0, 2, 0.0);

        nus_atm.Set_rel_error(1.0e-6);
        nus_atm.Set_abs_error(1.0e-6);
        nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);

        marray<double, 4> ini_state{nus_atm.GetNumCos(), nus_atm.GetNumE(), 2, 3};
        std::fill(ini_state.begin(), ini_state.end(), 0.0);
        for (size_t ci = 0; ci < nus_atm.GetNumCos(); ci++) {
            for (size_t ei = 0; ei < nus_atm.GetNumE(); ei++) {
                for (int rho = 0; rho < 2; rho++) {
                    ini_state[ci][ei][rho][1] = 1.0;
                }
            }
        }

        nus_atm.Set_initial_state(ini_state, flavor);
        nus_atm.Set_IncludeOscillations(true);

        nus_atm.EvolveState();
    }

    auto end = Clock::now();
    result.total_time_ms = Duration(end - start).count();
    result.time_per_iter_ms = result.total_time_ms / iterations;

    return result;
}

BenchmarkResult benchmark_atmospheric_glashow(int iterations = 2) {
    BenchmarkResult result;
    result.name = "Atm isoscalar (int.+Glashow)";
    result.description = "nuSQUIDSAtm: 10 GeV - 10 PeV (50 E) x 20 zenith, NC+CC+Glashow";
    result.iterations = iterations;
    result.passed = true;

    squids::Const units;

    double Emin = 10.0 * units.GeV;
    double Emax = 1.0e7 * units.GeV;  // 10 PeV to include Glashow resonance
    double czmin = -1.0;
    double czmax = 0.0;

    auto start = Clock::now();

    for (int iter = 0; iter < iterations; iter++) {
        nuSQUIDSAtm<> nus_atm(linspace(czmin, czmax, 20), logspace(Emin, Emax, 50), 3, both, true);

        nus_atm.Set_MixingAngle(0, 1, 0.563942);
        nus_atm.Set_MixingAngle(0, 2, 0.154085);
        nus_atm.Set_MixingAngle(1, 2, 0.785398);
        nus_atm.Set_SquareMassDifference(1, 7.65e-05);
        nus_atm.Set_SquareMassDifference(2, 0.00247);
        nus_atm.Set_CPPhase(0, 2, 0.0);

        nus_atm.Set_rel_error(1.0e-6);
        nus_atm.Set_abs_error(1.0e-6);
        nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);

        // Enable Glashow resonance
        nus_atm.Set_GlashowResonance(true);

        marray<double, 4> ini_state{nus_atm.GetNumCos(), nus_atm.GetNumE(), 2, 3};
        std::fill(ini_state.begin(), ini_state.end(), 0.0);
        for (size_t ci = 0; ci < nus_atm.GetNumCos(); ci++) {
            for (size_t ei = 0; ei < nus_atm.GetNumE(); ei++) {
                for (int rho = 0; rho < 2; rho++) {
                    ini_state[ci][ei][rho][1] = 1.0;
                }
            }
        }

        nus_atm.Set_initial_state(ini_state, flavor);
        nus_atm.Set_IncludeOscillations(true);

        nus_atm.EvolveState();
    }

    auto end = Clock::now();
    result.total_time_ms = Duration(end - start).count();
    result.time_per_iter_ms = result.total_time_ms / iterations;

    return result;
}

BenchmarkResult benchmark_atmospheric_full_physics(int iterations = 2) {
    BenchmarkResult result;
    result.name = "Atm isoscalar (full physics)";
    result.description = "nuSQUIDSAtm: 10 GeV - 10 PeV (50 E) x 20 zenith, NC+CC+Glashow+TauRegen";
    result.iterations = iterations;
    result.passed = true;

    squids::Const units;

    double Emin = 10.0 * units.GeV;
    double Emax = 1.0e7 * units.GeV;
    double czmin = -1.0;
    double czmax = 0.0;

    auto start = Clock::now();

    for (int iter = 0; iter < iterations; iter++) {
        nuSQUIDSAtm<> nus_atm(linspace(czmin, czmax, 20), logspace(Emin, Emax, 50), 3, both, true);

        nus_atm.Set_MixingAngle(0, 1, 0.563942);
        nus_atm.Set_MixingAngle(0, 2, 0.154085);
        nus_atm.Set_MixingAngle(1, 2, 0.785398);
        nus_atm.Set_SquareMassDifference(1, 7.65e-05);
        nus_atm.Set_SquareMassDifference(2, 0.00247);
        nus_atm.Set_CPPhase(0, 2, 0.0);

        nus_atm.Set_rel_error(1.0e-6);
        nus_atm.Set_abs_error(1.0e-6);
        nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);

        // Enable all physics
        nus_atm.Set_GlashowResonance(true);
        nus_atm.Set_TauRegeneration(true);

        marray<double, 4> ini_state{nus_atm.GetNumCos(), nus_atm.GetNumE(), 2, 3};
        std::fill(ini_state.begin(), ini_state.end(), 0.0);
        for (size_t ci = 0; ci < nus_atm.GetNumCos(); ci++) {
            for (size_t ei = 0; ei < nus_atm.GetNumE(); ei++) {
                for (int rho = 0; rho < 2; rho++) {
                    ini_state[ci][ei][rho][1] = 1.0;
                }
            }
        }

        nus_atm.Set_initial_state(ini_state, flavor);
        nus_atm.Set_IncludeOscillations(true);

        nus_atm.EvolveState();
    }

    auto end = Clock::now();
    result.total_time_ms = Duration(end - start).count();
    result.time_per_iter_ms = result.total_time_ms / iterations;

    return result;
}

BenchmarkResult benchmark_atmospheric_with_composition(int iterations = 2) {
    BenchmarkResult result;
    result.name = "Atm nuclear (PREM composition)";
    result.description = "nuSQUIDSAtm: 10 GeV - 1 PeV (50 E) x 20 zenith, PREM+composition";
    result.iterations = iterations;
    result.passed = true;

    squids::Const units;

    double Emin = 10.0 * units.GeV;
    double Emax = 1.0e6 * units.GeV;
    double czmin = -1.0;
    double czmax = 0.0;

    // Load Earth model with composition
    std::string prem_path = getResourcePath() + "/astro/EARTH_MODEL_PREM_wIso.dat";

    auto start = Clock::now();

    for (int iter = 0; iter < iterations; iter++) {
        nuSQUIDSAtm<> nus_atm(linspace(czmin, czmax, 20), logspace(Emin, Emax, 50), 3, both, true);

        // Use Earth model with nuclear composition
        auto earth = std::make_shared<EarthAtm>(prem_path);
        nus_atm.Set_EarthModel(earth);

        nus_atm.Set_MixingAngle(0, 1, 0.563942);
        nus_atm.Set_MixingAngle(0, 2, 0.154085);
        nus_atm.Set_MixingAngle(1, 2, 0.785398);
        nus_atm.Set_SquareMassDifference(1, 7.65e-05);
        nus_atm.Set_SquareMassDifference(2, 0.00247);
        nus_atm.Set_CPPhase(0, 2, 0.0);

        nus_atm.Set_rel_error(1.0e-6);
        nus_atm.Set_abs_error(1.0e-6);
        nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);

        marray<double, 4> ini_state{nus_atm.GetNumCos(), nus_atm.GetNumE(), 2, 3};
        std::fill(ini_state.begin(), ini_state.end(), 0.0);
        for (size_t ci = 0; ci < nus_atm.GetNumCos(); ci++) {
            for (size_t ei = 0; ei < nus_atm.GetNumE(); ei++) {
                for (int rho = 0; rho < 2; rho++) {
                    ini_state[ci][ei][rho][1] = 1.0;
                }
            }
        }

        nus_atm.Set_initial_state(ini_state, flavor);
        nus_atm.Set_IncludeOscillations(true);

        nus_atm.EvolveState();

        // Verify composition is being used by checking at a sample point
        auto& nus = nus_atm.GetnuSQuIDS(10);  // Get a middle zenith bin
        auto comp = earth->composition(*nus.GetTrack());
        if (comp.empty()) {
            result.passed = false;  // Composition should be non-empty
        }
    }

    auto end = Clock::now();
    result.total_time_ms = Duration(end - start).count();
    result.time_per_iter_ms = result.total_time_ms / iterations;

    return result;
}

int main(int argc, char* argv[]) {
    // Parse command line arguments for quick mode
    bool quick_mode = false;
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--quick" || std::string(argv[i]) == "-q") {
            quick_mode = true;
        }
    }

    print_header();

    std::vector<BenchmarkResult> results;

    // Run benchmarks - isoscalar mode (default, uses proton+neutron average)
    print_test_header("Single Energy Mode (isoscalar)");
    auto r1 = benchmark_single_energy(quick_mode ? 10 : 100);
    print_test_result(r1);
    results.push_back(r1);

    print_test_header("Multiple Energy Mode - isoscalar (no interactions)");
    auto r2 = benchmark_multiple_energy_no_interactions(quick_mode ? 2 : 10);
    print_test_result(r2);
    results.push_back(r2);

    print_test_header("Multiple Energy Mode - isoscalar (with interactions)");
    auto r3 = benchmark_multiple_energy_with_interactions(quick_mode ? 1 : 5);
    print_test_result(r3);
    results.push_back(r3);

    print_test_header("Atmospheric Mode - isoscalar (no interactions)");
    auto r4 = benchmark_atmospheric_no_interactions(quick_mode ? 1 : 3);
    print_test_result(r4);
    results.push_back(r4);

    print_test_header("Atmospheric Mode - isoscalar (with interactions)");
    auto r5 = benchmark_atmospheric_with_interactions(quick_mode ? 1 : 2);
    print_test_result(r5);
    results.push_back(r5);

    print_test_header("Atmospheric Mode - isoscalar (interactions + Glashow)");
    auto r6 = benchmark_atmospheric_glashow(quick_mode ? 1 : 2);
    print_test_result(r6);
    results.push_back(r6);

    print_test_header("Atmospheric Mode - isoscalar (full physics: int. + Glashow + Tau Regen)");
    auto r7 = benchmark_atmospheric_full_physics(quick_mode ? 1 : 2);
    print_test_result(r7);
    results.push_back(r7);

    // Nuclear composition mode (uses per-element cross sections with PREM composition)
    print_test_header("Atmospheric Mode - nuclear (PREM composition)");
    auto r8 = benchmark_atmospheric_with_composition(quick_mode ? 1 : 2);
    print_test_result(r8);
    results.push_back(r8);

    print_summary(results);

    return 0;
}
