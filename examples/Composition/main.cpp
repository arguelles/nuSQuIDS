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
 * This example demonstrates the use of body composition for nuclear cross sections.
 *
 * Bodies can provide composition information (elemental fractions) at each point
 * along a neutrino trajectory. This is useful for:
 * - Accurate cross section calculations with per-nucleus cross sections
 * - Studies of neutrino interactions in heterogeneous media
 * - Precise Earth matter effects with detailed PREM composition
 *
 * This example:
 * 1. Loads the PREM model with full nuclear composition
 * 2. Demonstrates accessing composition at various Earth depths
 * 3. Shows how to create a ConstantDensity body with custom composition
 */

#include <iostream>
#include <iomanip>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/resources.h>

using namespace nusquids;

// Helper function to get element name from PDGCode
std::string getElementName(PDGCode code) {
    switch(code) {
        case hydrogen:   return "H";
        case oxygen:     return "O";
        case sodium:     return "Na";
        case magnesium:  return "Mg";
        case aluminum:   return "Al";
        case silicon:    return "Si";
        case sulfur:     return "S";
        case calcium:    return "Ca";
        case iron:       return "Fe";
        case nickel:     return "Ni";
        default:         return "Unknown";
    }
}

int main()
{
    squids::Const units;

    std::cout << "========================================\n";
    std::cout << "  nuSQuIDS Body Composition Example\n";
    std::cout << "========================================\n\n";

    // =========================================================================
    // Part 1: Load Earth model with detailed nuclear composition
    // =========================================================================
    std::cout << "--- Part 1: Earth Model with Composition ---\n\n";

    // There are two PREM data files available:
    //
    // 1. EARTH_MODEL_PREM.dat (3 columns: r/R, density, ye)
    //    - Standard isoscalar model
    //    - Used with: EarthAtm() or EarthAtm("EARTH_MODEL_PREM.dat")
    //    - composition() returns empty map
    //
    // 2. EARTH_MODEL_PREM_wIso.dat (13 columns: r/R, density, ye, H, O, Na, Mg, Al, Si, S, Ca, Fe, Ni)
    //    - Model with full nuclear composition
    //    - Used with: EarthAtm("EARTH_MODEL_PREM_wIso.dat")
    //    - composition() returns element fractions at each depth
    //
    // The constructor automatically detects the file format based on column count.

    // Load the PREM model with isotopic/nuclear composition
    std::string prem_path = getResourcePath() + "/astro/EARTH_MODEL_PREM_wIso.dat";
    std::cout << "Loading Earth model WITH composition from: " << prem_path << "\n\n";

    auto earth = std::make_shared<EarthAtm>(prem_path);

    // For comparison, here's how to create an isoscalar Earth (no composition):
    // auto earth_isoscalar = std::make_shared<EarthAtm>();  // Uses default PREM
    // or: auto earth_isoscalar = std::make_shared<EarthAtm>(getResourcePath() + "/astro/EARTH_MODEL_PREM.dat");

    // Create tracks at different zenith angles to probe different depths
    std::vector<double> cos_zeniths = {-1.0, -0.8, -0.5, -0.2, 0.0};

    std::cout << "Composition at various depths (sampled at track midpoint):\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::setw(10) << "cos(zen)"
              << std::setw(12) << "depth(km)"
              << std::setw(12) << "density"
              << std::setw(8) << "Fe"
              << std::setw(8) << "O"
              << std::setw(8) << "Si"
              << std::setw(8) << "Mg"
              << "\n";
    std::cout << std::string(80, '-') << "\n";

    for (double cz : cos_zeniths) {
        auto track = earth->MakeTrackWithCosine(cz);

        // Sample at the midpoint of the track
        double x_mid = 0.5 * track.GetFinalX();
        track.SetX(x_mid);

        // Get density and composition at this point
        double density = earth->density(track);
        auto comp = earth->composition(track);

        // Calculate approximate depth (very rough)
        double depth_km = (earth->GetRadius() - x_mid / units.km);

        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << cz
                  << std::setw(12) << depth_km
                  << std::setw(12) << density;

        // Print some key element fractions
        std::cout << std::setprecision(4);
        std::cout << std::setw(8) << (comp.count(iron) ? comp.at(iron) : 0.0)
                  << std::setw(8) << (comp.count(oxygen) ? comp.at(oxygen) : 0.0)
                  << std::setw(8) << (comp.count(silicon) ? comp.at(silicon) : 0.0)
                  << std::setw(8) << (comp.count(magnesium) ? comp.at(magnesium) : 0.0)
                  << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // Part 2: Create ConstantDensity body with custom composition
    // =========================================================================
    std::cout << "--- Part 2: ConstantDensity with Custom Composition ---\n\n";

    // Create a water-equivalent material (useful for ice/water detectors)
    std::map<PDGCode, double> water_composition;
    water_composition[hydrogen] = 2.0/3.0;  // 2 H atoms per water molecule
    water_composition[oxygen] = 1.0/3.0;    // 1 O atom per water molecule

    double water_density = 1.0;  // g/cm^3
    double water_ye = 0.555;     // electron fraction for water

    auto water = std::make_shared<ConstantDensity>(water_density, water_ye, water_composition);
    auto water_track = std::make_shared<ConstantDensity::Track>(1000.0 * units.km);

    std::cout << "Water-equivalent material:\n";
    std::cout << "  Density: " << water->density(*water_track) << " g/cm^3\n";
    std::cout << "  Electron fraction: " << water->ye(*water_track) << "\n";
    std::cout << "  Composition:\n";
    auto water_comp = water->composition(*water_track);
    for (const auto& elem : water_comp) {
        std::cout << "    " << getElementName(elem.first) << ": "
                  << std::fixed << std::setprecision(4) << elem.second << "\n";
    }
    std::cout << "\n";

    // =========================================================================
    // Part 3: Run a simple propagation with composition-aware body
    // =========================================================================
    std::cout << "--- Part 3: Propagation with Composition ---\n\n";

    // Create nuSQuIDS instance
    nuSQUIDS nus(logspace(1.0*units.GeV, 1.0e6*units.GeV, 100), 3, neutrino, true);

    // Set mixing parameters
    nus.Set_MixingAngle(0,1, 0.5837);  // theta_12
    nus.Set_MixingAngle(0,2, 0.1496);  // theta_13
    nus.Set_MixingAngle(1,2, 0.8552);  // theta_23
    nus.Set_SquareMassDifference(1, 7.42e-5);  // dm21^2
    nus.Set_SquareMassDifference(2, 2.51e-3);  // dm31^2

    // Use Earth model with composition - must set body/track BEFORE initial state
    nus.Set_Body(earth);
    auto prop_track = std::make_shared<EarthAtm::Track>(earth->MakeTrackWithCosine(-1.0));  // Vertical upgoing
    nus.Set_Track(prop_track);

    // Set initial state (power law flux)
    marray<double, 2> inistate({100, 3});
    for (size_t ie = 0; ie < 100; ie++) {
        for (size_t iflv = 0; iflv < 3; iflv++) {
            inistate[ie][iflv] = (iflv == 1) ? 1.0 : 0.0;  // Pure muon neutrino
        }
    }
    nus.Set_initial_state(inistate, flavor);

    std::cout << "Propagating muon neutrinos through Earth (cos_zen = -1)...\n";
    nus.EvolveState();
    std::cout << "Done!\n\n";

    // Print survival probability at a few energies
    std::cout << "Muon neutrino survival probability:\n";
    std::cout << std::setw(15) << "Energy (GeV)" << std::setw(15) << "P(nu_mu->nu_mu)\n";
    auto e_range = nus.GetERange();
    std::vector<double> sample_energies = {10.0, 100.0, 1000.0, 10000.0};
    for (double E : sample_energies) {
        // Find energy index
        size_t ie = 0;
        for (size_t i = 0; i < e_range.size(); i++) {
            if (e_range[i] >= E * units.GeV) {
                ie = i;
                break;
            }
        }
        double prob = nus.EvalFlavorAtNode(1, ie);  // nu_mu survival (flavor=1)
        std::cout << std::setw(15) << E << std::setw(15) << std::setprecision(4) << prob << "\n";
    }

    std::cout << "\n========================================\n";
    std::cout << "  Example Complete!\n";
    std::cout << "========================================\n";

    return 0;
}
