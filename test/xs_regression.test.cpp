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

/*******************************************************************************
 * Cross Section Regression Test
 *
 * This test verifies that neutrino propagation results remain consistent
 * across code changes for:
 *
 * 1. Isoscalar cross sections (single target type)
 * 2. Proton/neutron cross sections (two target types)
 * 3. Nuclear cross sections (multiple nuclear targets with composition)
 *
 * Physics configurations tested:
 * - Oscillations + Interactions (no Glashow)
 * - Oscillations + Interactions + Glashow (except nuclear case baseline from nuclearxs-dev)
 *
 * The test uses an E^-2 power law spectrum propagated through Earth.
 ******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/xsections.h>
#include <nuSQuIDS/resources.h>

using namespace nusquids;

// Test configuration
const unsigned int numneu = 3;
const unsigned int n_energies = 40;  // 40 nodes for faster execution
const double E_min = 1.0e2;   // 100 GeV
const double E_max = 1.0e7;   // 10 PeV (5 decades)

// Print flux values at selected energies for comparison
void print_fluxes(const std::string& label, nuSQUIDS& nus, const squids::Const& units) {
    auto e_range = nus.GetERange();

    // Sample at 5 points across the energy range
    std::vector<size_t> sample_indices = {0, 10, 20, 30, 39};

    std::cout << std::scientific << std::setprecision(6);
    for (size_t idx : sample_indices) {
        if (idx >= e_range.size()) continue;
        double E = e_range[idx];

        std::cout << label << " E=" << E/units.GeV << "GeV";
        for (unsigned int flv = 0; flv < numneu; flv++) {
            // neutrino
            double flux_nu = nus.EvalFlavorAtNode(flv, idx, 0);
            std::cout << " nu" << flv << "=" << flux_nu;
        }
        for (unsigned int flv = 0; flv < numneu; flv++) {
            // antineutrino
            double flux_nubar = nus.EvalFlavorAtNode(flv, idx, 1);
            std::cout << " nubar" << flv << "=" << flux_nubar;
        }
        std::cout << std::endl;
    }
}

// Run propagation test with given cross sections and physics options
void run_test(const std::string& label,
              CrossSectionLibrary xs,
              std::shared_ptr<Body> body,
              std::shared_ptr<Body::Track> track,
              bool enable_glashow) {

    squids::Const units;

    // Create energy array
    auto energies = logspace(E_min * units.GeV, E_max * units.GeV, n_energies);

    // Create nuSQuIDS instance - wrap in shared_ptr
    auto xs_ptr = std::make_shared<CrossSectionLibrary>(std::move(xs));
    nuSQUIDS nus(energies, numneu, both, true, xs_ptr);

    // Set mixing parameters (standard values)
    nus.Set_MixingAngle(0, 1, 0.5837);   // theta_12
    nus.Set_MixingAngle(0, 2, 0.1496);   // theta_13
    nus.Set_MixingAngle(1, 2, 0.8552);   // theta_23
    nus.Set_SquareMassDifference(1, 7.42e-5);  // dm21^2
    nus.Set_SquareMassDifference(2, 2.51e-3);  // dm31^2
    nus.Set_CPPhase(0, 2, 0.0);

    // Set physics options
    nus.Set_GlashowResonance(enable_glashow);
    nus.Set_TauRegeneration(false);

    // Set integration parameters for reproducibility
    nus.Set_rel_error(1.0e-10);
    nus.Set_abs_error(1.0e-10);

    // Set body and track
    nus.Set_Body(body);
    nus.Set_Track(track);

    // Set initial state: E^-2 power law for all flavors (both nu and nubar)
    marray<double, 3> inistate({n_energies, 2, numneu});  // [energy, rho, flavor]
    for (size_t ie = 0; ie < n_energies; ie++) {
        double E = energies[ie];
        double flux = pow(E / units.GeV, -2.0);  // E^-2 spectrum
        for (size_t rho = 0; rho < 2; rho++) {   // 0=nu, 1=nubar
            for (size_t flv = 0; flv < numneu; flv++) {
                inistate[ie][rho][flv] = flux;
            }
        }
    }
    nus.Set_initial_state(inistate, flavor);

    // Evolve
    nus.EvolveState();

    // Print results
    print_fluxes(label, nus, units);
}

int main() {
    squids::Const units;

    std::cout << "=== Cross Section Regression Test ===" << std::endl;
    std::cout << std::endl;

    // Create Earth body with standard PREM (no composition)
    auto earth_std = std::make_shared<EarthAtm>();
    auto track_std = std::make_shared<EarthAtm::Track>(earth_std->MakeTrackWithCosine(-0.5));

    // Create Earth body with composition for nuclear XS
    std::string prem_comp_path = getResourcePath() + "/astro/EARTH_MODEL_PREM_wIso.dat";
    auto earth_comp = std::make_shared<EarthAtm>(prem_comp_path);
    auto track_comp = std::make_shared<EarthAtm::Track>(earth_comp->MakeTrackWithCosine(-0.5));

    // =========================================================================
    // Test 1: Isoscalar cross sections
    // =========================================================================
    std::cout << "--- Isoscalar Cross Sections ---" << std::endl;
    {
        std::string xsdir = getResourcePath() + "/xsections/";

        // Create isoscalar cross section library
        CrossSectionLibrary xs_isoscalar;
        xs_isoscalar.addTarget(isoscalar_nucleon,
            NeutrinoDISCrossSectionsFromTables(xsdir + "csms_square.h5"));
        xs_isoscalar.addTarget(electron, GlashowResonanceCrossSection());

        // Without Glashow
        run_test("ISOSCALAR_NOGR", xs_isoscalar, earth_std, track_std, false);

        // Create again for second test (since we move it)
        CrossSectionLibrary xs_isoscalar2;
        xs_isoscalar2.addTarget(isoscalar_nucleon,
            NeutrinoDISCrossSectionsFromTables(xsdir + "csms_square.h5"));
        xs_isoscalar2.addTarget(electron, GlashowResonanceCrossSection());

        // With Glashow
        run_test("ISOSCALAR_GR", xs_isoscalar2, earth_std, track_std, true);
    }
    std::cout << std::endl;

    // =========================================================================
    // Test 2: Proton/Neutron cross sections (default)
    // =========================================================================
    std::cout << "--- Proton/Neutron Cross Sections ---" << std::endl;
    {
        // Without Glashow
        run_test("PN_NOGR", loadDefaultCrossSections(), earth_std, track_std, false);

        // With Glashow
        run_test("PN_GR", loadDefaultCrossSections(), earth_std, track_std, true);
    }
    std::cout << std::endl;

    // =========================================================================
    // Test 3: Nuclear cross sections with composition
    // =========================================================================
    std::cout << "--- Nuclear Cross Sections ---" << std::endl;
    {
        // Without Glashow (baseline from nuclearxs-dev)
        run_test("NUCLEAR_NOGR", loadWCG24NuclearCrossSections(), earth_comp, track_comp, false);

        // With Glashow (new feature, no baseline from nuclearxs-dev)
        // We include this for future regression tracking
        run_test("NUCLEAR_GR", loadWCG24NuclearCrossSections(), earth_comp, track_comp, true);
    }
    std::cout << std::endl;

    std::cout << "=== Test Complete ===" << std::endl;

    return 0;
}
