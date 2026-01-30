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
*      Connor Sponsler                                                        *
*      Carlos Arguelles (University of Wisconsin Madison)                     *
*         carguelles@icecube.wisc.edu                                         *
******************************************************************************/

/*
 * This example demonstrates how to use EmittingEarthAtm to simulate
 * atmospheric neutrino production and propagation through the Earth.
 *
 * EmittingEarthAtm injects neutrino flux during propagation based on
 * a production profile (e.g., from MCEq). This allows simulating the
 * continuous production of atmospheric neutrinos as cosmic ray showers
 * develop in the atmosphere.
 *
 * Before running this example, generate the production profile using:
 *   cd data/atmos_prod
 *   python generate_production_profile.py
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <nuSQuIDS/nuSQuIDS.h>

using namespace nusquids;

int main()
{
  squids::Const units;

  // Set up energy range (1 GeV to 10 TeV)
  unsigned int n_energies = 100;
  marray<double,1> E_range = logspace(1.0*units.GeV, 10000.0*units.GeV, n_energies);

  // Set up zenith angle (upgoing neutrinos, cos(zenith) = -1)
  double coszen = -0.7;  // Neutrinos coming from below horizon

  // Create EmittingEarthAtm body
  // This will load the default production profile from data/atmos_prod/PROD_MODEL_MCEQ.dat
  std::shared_ptr<EmittingEarthAtm> emitting_earth = std::make_shared<EmittingEarthAtm>();

  // Create track through the atmosphere
  auto track = std::make_shared<EarthAtm::Track>(emitting_earth->MakeTrackWithCosine(coszen));

  // Create nuSQUIDS object with both neutrinos and antineutrinos
  nuSQUIDS nus(E_range, 3, both, true);

  // Set the emitting Earth body and track
  nus.Set_Body(emitting_earth);
  nus.Set_Track(track);

  // Set oscillation parameters
  nus.Set_MixingAngle(0, 1, 0.5836);       // theta_12
  nus.Set_MixingAngle(0, 2, 0.1496);       // theta_13
  nus.Set_MixingAngle(1, 2, 0.7854);       // theta_23
  nus.Set_SquareMassDifference(1, 7.5e-5); // Delta m^2_21 [eV^2]
  nus.Set_SquareMassDifference(2, 2.5e-3); // Delta m^2_31 [eV^2]
  nus.Set_CPPhase(0, 2, 0.0);              // delta_CP

  // Enable neutrino sources (this is crucial for EmittingEarthAtm to work!)
  nus.Set_NeutrinoSources(true);

  // Set numerical precision
  nus.Set_rel_error(1.0e-6);
  nus.Set_abs_error(1.0e-6);
  nus.Set_GSL_step(gsl_odeiv2_step_rk4);

  // Show progress
  nus.Set_ProgressBar(true);

  // Initialize with zero flux (neutrinos will be injected during propagation)
  marray<double,3> inistate{n_energies, 2, 3};  // [energy][rho][flavor]
  std::fill(inistate.begin(), inistate.end(), 0.0);
  nus.Set_initial_state(inistate, flavor);

  std::cout << "Propagating atmospheric neutrinos with continuous injection..." << std::endl;
  std::cout << "  cos(zenith) = " << coszen << std::endl;
  std::cout << "  Energy range: " << E_range[0]/units.GeV << " - " << E_range[n_energies-1]/units.GeV << " GeV" << std::endl;

  // Evolve the system
  nus.EvolveState();

  // Write output
  std::cout << "\nWriting fluxes to fluxes_flavor.txt..." << std::endl;
  std::ofstream file("fluxes_flavor.txt");
  file << "# E[GeV] nue_flux numu_flux nutau_flux nuebar_flux numubar_flux nutaubar_flux" << std::endl;

  for(unsigned int ie = 0; ie < n_energies; ie++){
    double E = E_range[ie];
    file << E/units.GeV;
    // Neutrino flavors
    for(unsigned int fl = 0; fl < 3; fl++){
      file << " " << nus.EvalFlavor(fl, E, 0);  // rho=0 is neutrino
    }
    // Antineutrino flavors
    for(unsigned int fl = 0; fl < 3; fl++){
      file << " " << nus.EvalFlavor(fl, E, 1);  // rho=1 is antineutrino
    }
    file << std::endl;
  }
  file.close();

  std::cout << "Done!" << std::endl;
  return 0;
}
