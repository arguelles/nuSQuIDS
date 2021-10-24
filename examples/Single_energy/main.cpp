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

#include <vector>
#include <iostream>
#include "nuSQuIDS/nuSQuIDS.h"
/*
 * This file demonstrates how to calculate neutrino oscillation
 * probabilities for different initial flavor states, mixing angles
 * and bodies in the single energy mode.
 */

using namespace nusquids;

int main()
{
  // We must first create a nuSQUIDS object. In order to do this
  // we must specify the number of neutrino flavors and if we are
  // going to consider neutrino or antineutrino oscillations.

  // In this example we set N_neutrino = 3 and Type = "neutrino".
  nuSQUIDS nus(3,neutrino);

  // We will now set the mixing parameters and square mass
  // differences. The angles will be given in radiant and the
  // differences in eV^2.

  // We use the standard parametrization as described in the
  // documentation.

  // mixing angles
  nus.Set_MixingAngle(0,1,0.563942);
  nus.Set_MixingAngle(0,2,0.154085);
  nus.Set_MixingAngle(1,2,0.785398);
  // square mass differences
  nus.Set_SquareMassDifference(1,7.65e-05);
  nus.Set_SquareMassDifference(2,0.00247);
  // CP phase
  nus.Set_CPPhase(0,2,0.0);

  // Define a Const object for handleing the units
  squids::Const units;

  // Now we set the neutrino energy which we are interested on.
  // Energies are always given in natural units. To handle
  // the units the nuSQUIDS object has a unit subclass
  // which contains the most used units.

  nus.Set_E(10.0*units.GeV);

  // To calculate atmospheric neutrino oscillation probabilities
  // we need to specify a different body.

  std::cout << "*** Earth Atmospheric Neutrino Osc ***" << std::endl;
  //here we set the zenith angle the object ant the track, that is basically parametrized with the zenith angle
  double phi = acos(-1.0);
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  auto earth_atm_track = std::make_shared<EarthAtm::Track>(earth_atm->MakeTrack(phi));  
  nus.Set_Body(earth_atm);
  nus.Set_Track(earth_atm_track);

  //Here we set the initial state for the flavor, a muon state for this case
  marray<double,1> ini_state({3},{0,1,0});
  nus.Set_initial_state(ini_state,flavor);

  
  // setup integration settings
  nus.Set_rel_error(1.0e-20);
  nus.Set_abs_error(1.0e-20);
  
  // We can change the energy
  nus.Set_E(100.0*units.GeV);
  
  // We reset the initial condition
  nus.Set_initial_state(ini_state,flavor);
  
  std::cout << "In state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }
  
  // We do the calculation
  nus.EvolveState();
  
  // Output the result
  std::cout << "Out state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
