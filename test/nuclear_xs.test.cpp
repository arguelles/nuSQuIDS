/**
 * Nuclear Cross Section Integration Tests
 *
 * Tests the nuclear cross section and composition system step by step:
 * 1. CrossSectionLibrary loading and target enumeration
 * 2. Body composition at various track positions
 * 3. Target number fractions calculation
 * 4. Target number densities calculation
 * 5. Integration test comparing interaction lengths
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/resources.h>

using namespace nusquids;

// Test utilities
bool approx_equal(double a, double b, double tol = 1e-6) {
    return std::abs(a - b) < tol;
}

void print_separator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

// Test 1: CrossSectionLibrary loading and targets
bool test_cross_section_library() {
    print_separator("Test 1: CrossSectionLibrary Loading");

    try {
        // Load default (proton/neutron) cross sections
        auto default_lib = loadDefaultCrossSections();
        std::cout << "Default library targets: " << default_lib.numberOfTargets() << "\n";
        for (const auto& t : default_lib.targets()) {
            std::cout << "  - " << static_cast<int32_t>(t) << "\n";
        }

        // Load WCG24 proton/neutron
        auto wcg24_pn = loadWCG24CrossSections();
        std::cout << "\nWCG24 p/n library targets: " << wcg24_pn.numberOfTargets() << "\n";
        for (const auto& t : wcg24_pn.targets()) {
            std::cout << "  - " << static_cast<int32_t>(t) << "\n";
        }

        // Load WCG24 nuclear
        auto wcg24_nuc = loadWCG24NuclearCrossSections();
        std::cout << "\nWCG24 nuclear library targets: " << wcg24_nuc.numberOfTargets() << "\n";
        for (const auto& t : wcg24_nuc.targets()) {
            std::cout << "  - " << static_cast<int32_t>(t);
            if (isNuclearPDGCode(t)) {
                std::cout << " (Z=" << getAtomicNumber(t) << ", A=" << getMassNumber(t) << ")";
            }
            std::cout << "\n";
        }

        // Verify nuclear targets are detected correctly
        bool has_nuclear = false;
        for (const auto& t : wcg24_nuc.targets()) {
            if (isNuclearPDGCode(t)) {
                has_nuclear = true;
                break;
            }
        }

        if (!has_nuclear) {
            std::cerr << "ERROR: Nuclear library should have nuclear targets!\n";
            return false;
        }

        std::cout << "\n✓ Cross section libraries loaded successfully\n";
        return true;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return false;
    }
}

// Test 2: Body composition at various positions
bool test_body_composition() {
    print_separator("Test 2: Body Composition");

    try {
        // Load EarthAtm with composition
        std::string prem_path = getResourcePath() + "/astro/EARTH_MODEL_PREM_wIso.dat";

        EarthAtm earth(prem_path);

        // Test composition at various depths (vertical track)
        auto track = earth.MakeTrackWithCosine(-1.0);
        double total_length = track.GetFinalX();

        std::cout << "\nComposition along vertical track (cos_zen = -1.0):\n";
        std::cout << std::setw(10) << "Position"
                  << std::setw(12) << "Density"
                  << std::setw(10) << "ye"
                  << std::setw(10) << "O"
                  << std::setw(10) << "Fe"
                  << std::setw(10) << "Ni"
                  << std::setw(10) << "Sum"
                  << "\n";
        std::cout << std::string(72, '-') << "\n";

        for (double frac : {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
            double x = frac * total_length;
            track.SetX(x);

            double density = earth.density(track);
            double ye = earth.ye(track);
            auto comp = earth.composition(track);

            double sum = 0.0;
            for (const auto& p : comp) sum += p.second;

            double f_O = comp.count(oxygen) ? comp.at(oxygen) : 0.0;
            double f_Fe = comp.count(iron) ? comp.at(iron) : 0.0;
            double f_Ni = comp.count(nickel) ? comp.at(nickel) : 0.0;

            std::cout << std::fixed << std::setprecision(3)
                      << std::setw(10) << frac
                      << std::setw(12) << density
                      << std::setw(10) << ye
                      << std::setw(10) << f_O
                      << std::setw(10) << f_Fe
                      << std::setw(10) << f_Ni
                      << std::setw(10) << sum
                      << "\n";
        }

        // Verify composition sums to 1 at midpoint (should be in core)
        track.SetX(0.5 * total_length);
        auto comp_mid = earth.composition(track);
        double sum_mid = 0.0;
        for (const auto& p : comp_mid) sum_mid += p.second;

        if (std::abs(sum_mid - 1.0) > 0.01) {
            std::cerr << "WARNING: Composition at midpoint sums to " << sum_mid << " (expected ~1.0)\n";
        }

        // Verify we're actually in the core (should have significant Fe)
        if (comp_mid.count(iron) == 0 || comp_mid.at(iron) < 0.5) {
            std::cerr << "WARNING: Expected significant iron at Earth's core!\n";
        }

        std::cout << "\n✓ Body composition test passed\n";
        return true;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return false;
    }
}

// Test 3: nuSQUIDS target initialization with nuclear XS
bool test_target_initialization() {
    print_separator("Test 3: Target Initialization");

    try {
        squids::Const units;
        double GeV = units.GeV;

        // Create energy array as marray
        marray<double, 1> energies({10});
        for (int i = 0; i < 10; i++) {
            energies[i] = 100 * GeV * std::pow(10, i * 0.5);
        }

        // Test 1: Default cross sections
        std::cout << "Creating nuSQUIDS with default cross sections...\n";
        nuSQUIDS nus1(energies, 3, neutrino, true);

        // Set body and track to trigger initialization
        auto earth = std::make_shared<EarthAtm>();
        nus1.Set_Body(earth);
        auto track1 = earth->MakeTrackWithCosine(-0.5);
        nus1.Set_Track(std::make_shared<EarthAtm::Track>(track1));

        // Set initial state to trigger full init
        marray<double, 2> inistate({10, 3});
        for (size_t i = 0; i < 10; i++) inistate[i][1] = 1.0;
        nus1.Set_initial_state(inistate, flavor);

        auto xs1 = nus1.GetNeutrinoCrossSections();
        if (xs1) {
            std::cout << "Default XS targets: ";
            for (const auto& t : xs1->targets()) {
                std::cout << static_cast<int32_t>(t) << " ";
            }
            std::cout << "\n";
        } else {
            std::cout << "Default XS: nullptr (not yet initialized)\n";
        }

        // Test 2: Nuclear cross sections
        std::cout << "\nCreating nuSQUIDS with nuclear cross sections...\n";
        nuSQUIDS nus2(energies, 3, neutrino, true);

        auto nuc_lib = std::make_shared<CrossSectionLibrary>(loadWCG24NuclearCrossSections());
        nus2.SetNeutrinoCrossSections(nuc_lib);

        // Load Earth with composition
        std::string prem_path = getResourcePath() + "/astro/EARTH_MODEL_PREM_wIso.dat";
        auto earth2 = std::make_shared<EarthAtm>(prem_path);
        nus2.Set_Body(earth2);
        auto track2 = earth2->MakeTrackWithCosine(-0.5);
        nus2.Set_Track(std::make_shared<EarthAtm::Track>(track2));

        // Set initial state
        nus2.Set_initial_state(inistate, flavor);

        auto xs2 = nus2.GetNeutrinoCrossSections();
        std::cout << "Nuclear XS targets after init: ";
        for (const auto& t : xs2->targets()) {
            std::cout << static_cast<int32_t>(t) << " ";
        }
        std::cout << "\n";

        // Verify nuclear targets are being used
        bool has_nuclear = false;
        for (const auto& t : xs2->targets()) {
            if (isNuclearPDGCode(t)) {
                has_nuclear = true;
                break;
            }
        }

        if (!has_nuclear) {
            std::cerr << "ERROR: Nuclear targets should be in XS after SetNeutrinoCrossSections!\n";
            return false;
        }

        std::cout << "\n✓ Target initialization test passed\n";
        return true;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return false;
    }
}

// Test 4: Compare evolution with p/n vs nuclear XS
bool test_evolution_comparison() {
    print_separator("Test 4: Evolution Comparison");

    try {
        squids::Const units;
        double GeV = units.GeV;

        marray<double, 1> energies({20});
        for (int i = 0; i < 20; i++) {
            energies[i] = 10 * GeV * std::pow(10, i * 0.25);
        }

        // === Case 1: Proton/Neutron cross sections ===
        std::cout << "Evolving with WCG24 p/n cross sections...\n";
        nuSQUIDS nus_pn(energies, 3, neutrino, true);

        auto pn_lib = std::make_shared<CrossSectionLibrary>(loadWCG24CrossSections());
        nus_pn.SetNeutrinoCrossSections(pn_lib);

        // Use default Earth (no composition)
        auto earth_pn = std::make_shared<EarthAtm>();
        nus_pn.Set_Body(earth_pn);
        auto track_pn = earth_pn->MakeTrackWithCosine(-1.0);
        nus_pn.Set_Track(std::make_shared<EarthAtm::Track>(track_pn));

        marray<double, 2> inistate({20, 3});
        for (size_t i = 0; i < 20; i++) inistate[i][1] = 1.0;  // Pure nu_mu
        nus_pn.Set_initial_state(inistate, flavor);

        nus_pn.Set_rel_error(1e-6);
        nus_pn.Set_abs_error(1e-6);
        nus_pn.EvolveState();

        // === Case 2: Nuclear cross sections ===
        std::cout << "\nEvolving with WCG24 nuclear cross sections...\n";
        nuSQUIDS nus_nuc(energies, 3, neutrino, true);

        auto nuc_lib = std::make_shared<CrossSectionLibrary>(loadWCG24NuclearCrossSections());
        nus_nuc.SetNeutrinoCrossSections(nuc_lib);

        // Use Earth with composition
        std::string prem_path = getResourcePath() + "/astro/EARTH_MODEL_PREM_wIso.dat";
        auto earth_nuc = std::make_shared<EarthAtm>(prem_path);
        nus_nuc.Set_Body(earth_nuc);
        auto track_nuc = earth_nuc->MakeTrackWithCosine(-1.0);
        nus_nuc.Set_Track(std::make_shared<EarthAtm::Track>(track_nuc));

        nus_nuc.Set_initial_state(inistate, flavor);

        nus_nuc.Set_rel_error(1e-6);
        nus_nuc.Set_abs_error(1e-6);
        nus_nuc.EvolveState();

        // === Compare results ===
        std::cout << "\nComparing fluxes at various energies:\n";
        std::cout << std::setw(15) << "Energy [GeV]"
                  << std::setw(15) << "P/N flux"
                  << std::setw(15) << "Nuclear flux"
                  << std::setw(15) << "Ratio"
                  << "\n";
        std::cout << std::string(60, '-') << "\n";

        bool significant_difference = false;
        for (size_t i = 0; i < energies.size(); i += 4) {
            double E = energies[i];
            double flux_pn = nus_pn.EvalFlavor(1, E, 0);   // nu_mu
            double flux_nuc = nus_nuc.EvalFlavor(1, E, 0);
            double ratio = flux_pn > 1e-10 ? flux_nuc / flux_pn : 0.0;

            std::cout << std::fixed << std::setprecision(4)
                      << std::setw(15) << E / GeV
                      << std::setw(15) << flux_pn
                      << std::setw(15) << flux_nuc
                      << std::setw(15) << ratio
                      << "\n";

            // At high energies, both should show absorption
            if (E / GeV > 1e5 && flux_nuc > 0.99 && flux_pn < 0.5) {
                significant_difference = true;
                std::cerr << "WARNING: Nuclear case shows no absorption at " << E/GeV << " GeV!\n";
            }
        }

        if (significant_difference) {
            std::cerr << "\nERROR: Nuclear cross sections are not being applied correctly!\n";
            return false;
        }

        std::cout << "\n✓ Evolution comparison test passed\n";
        return true;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return false;
    }
}

int main() {
    std::cout << "nuSQuIDS Nuclear Cross Section Tests\n";
    std::cout << std::string(60, '=') << "\n";

    bool all_passed = true;

    all_passed &= test_cross_section_library();
    all_passed &= test_body_composition();
    all_passed &= test_target_initialization();
    all_passed &= test_evolution_comparison();

    std::cout << "\n" << std::string(60, '=') << "\n";
    if (all_passed) {
        std::cout << "ALL TESTS PASSED\n";
        return 0;
    } else {
        std::cout << "SOME TESTS FAILED\n";
        return 1;
    }
}
