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
 * EarthComposition Example - Comparing Isoscalar vs Nuclear Composition
 *
 * This example demonstrates the difference between:
 * 1. Isoscalar mode: Standard proton+neutron averaged cross sections
 * 2. Nuclear composition mode: Per-element cross sections with PREM composition
 *
 * We propagate atmospheric neutrinos through Earth using nuSQUIDSAtm and compare
 * the final flux ratios for both modes.
 *
 * Output: flux_ratios_isoscalar.dat and flux_ratios_composition.dat
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/resources.h>

using namespace nusquids;

void run_propagation(bool use_composition, const std::string& output_filename) {
    squids::Const units;

    // Energy range: 10 GeV to 1 PeV
    double Emin = 10.0 * units.GeV;
    double Emax = 1.0e6 * units.GeV;
    unsigned int n_energies = 100;

    // Zenith angles: upgoing through horizontal
    double cos_zenith_min = -1.0;
    double cos_zenith_max = 0.0;
    unsigned int n_cos_zenith = 20;

    auto cos_zeniths = linspace(cos_zenith_min, cos_zenith_max, n_cos_zenith);
    auto energies = logspace(Emin, Emax, n_energies);

    // Create nuSQUIDSAtm with interactions
    nuSQUIDSAtm<> nus_atm(cos_zeniths, energies, 3, both, true);

    // Set Earth model
    if (use_composition) {
        // Load PREM with nuclear composition (13-column file)
        std::string prem_path = getResourcePath() + "/astro/EARTH_MODEL_PREM_wIso.dat";
        auto earth = std::make_shared<EarthAtm>(prem_path);
        nus_atm.Set_EarthModel(earth);
        std::cout << "Using PREM with nuclear composition" << std::endl;
    } else {
        // Use default isoscalar PREM (3-column file)
        auto earth = std::make_shared<EarthAtm>();
        nus_atm.Set_EarthModel(earth);
        std::cout << "Using isoscalar PREM" << std::endl;
    }

    // Set mixing parameters (NuFit 5.0 values)
    nus_atm.Set_MixingAngle(0, 1, 0.5837);   // theta_12
    nus_atm.Set_MixingAngle(0, 2, 0.1496);   // theta_13
    nus_atm.Set_MixingAngle(1, 2, 0.8552);   // theta_23
    nus_atm.Set_SquareMassDifference(1, 7.42e-5);   // dm21^2
    nus_atm.Set_SquareMassDifference(2, 2.51e-3);   // dm31^2
    nus_atm.Set_CPPhase(0, 2, 0.0);

    // Set numerical parameters
    nus_atm.Set_rel_error(1.0e-6);
    nus_atm.Set_abs_error(1.0e-6);
    nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);

    // Set initial state: E^-2 spectrum (atmospheric-like)
    // Initial state shape: [n_cos_zenith, n_energies, 2 (nu/nubar), 3 (flavors)]
    marray<double, 4> initial_state{n_cos_zenith, n_energies, 2, 3};
    std::fill(initial_state.begin(), initial_state.end(), 0.0);

    for (unsigned int icz = 0; icz < n_cos_zenith; icz++) {
        for (unsigned int ie = 0; ie < n_energies; ie++) {
            double E = energies[ie] / units.GeV;
            double flux = std::pow(E, -2.0);  // E^-2 spectrum

            for (int rho = 0; rho < 2; rho++) {  // neutrino and antineutrino
                // Atmospheric ratio approximately: nu_e : nu_mu : nu_tau = 1 : 2 : 0
                initial_state[icz][ie][rho][0] = flux * 1.0;  // nu_e
                initial_state[icz][ie][rho][1] = flux * 2.0;  // nu_mu
                initial_state[icz][ie][rho][2] = flux * 0.0;  // nu_tau
            }
        }
    }

    nus_atm.Set_initial_state(initial_state, flavor);
    nus_atm.Set_IncludeOscillations(true);

    // Evolve
    std::cout << "Propagating..." << std::endl;
    nus_atm.EvolveState();
    std::cout << "Done!" << std::endl;

    // Write output file
    std::ofstream outfile(output_filename);
    outfile << "# EarthComposition Example Output" << std::endl;
    outfile << "# " << (use_composition ? "Nuclear composition mode" : "Isoscalar mode") << std::endl;
    outfile << "# Columns: cos_zenith  E_GeV  ratio_nue  ratio_numu  ratio_nutau" << std::endl;
    outfile << "# ratio = (final_flux) / (initial_flux), summed over nu+nubar" << std::endl;
    outfile << std::scientific << std::setprecision(6);

    for (unsigned int icz = 0; icz < n_cos_zenith; icz++) {
        double cz = cos_zeniths[icz];
        for (unsigned int ie = 0; ie < n_energies; ie++) {
            double E = energies[ie];
            double E_GeV = E / units.GeV;

            outfile << std::setw(14) << cz << " " << std::setw(14) << E_GeV;

            for (int flv = 0; flv < 3; flv++) {
                // Sum over neutrino and antineutrino
                double final_flux = nus_atm.EvalFlavor(flv, cz, E, 0) +
                                   nus_atm.EvalFlavor(flv, cz, E, 1);
                double initial_flux = initial_state[icz][ie][0][flv] +
                                     initial_state[icz][ie][1][flv];

                double ratio = (initial_flux > 0) ? final_flux / initial_flux : 0.0;
                outfile << " " << std::setw(14) << ratio;
            }
            outfile << std::endl;
        }
        outfile << std::endl;  // Blank line between zenith angles (for gnuplot)
    }

    outfile.close();
    std::cout << "Wrote: " << output_filename << std::endl;
}

int main()
{
    std::cout << "========================================" << std::endl;
    std::cout << "  EarthComposition Example" << std::endl;
    std::cout << "  Isoscalar vs Nuclear Composition" << std::endl;
    std::cout << "========================================" << std::endl << std::endl;

    // Run isoscalar propagation
    run_propagation(false, "flux_ratios_isoscalar.dat");

    std::cout << std::endl;

    // Run composition propagation
    run_propagation(true, "flux_ratios_composition.dat");

    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  Example Complete!" << std::endl;
    std::cout << "  Output files:" << std::endl;
    std::cout << "    - flux_ratios_isoscalar.dat" << std::endl;
    std::cout << "    - flux_ratios_composition.dat" << std::endl;
    std::cout << "========================================" << std::endl;

    return 0;
}
