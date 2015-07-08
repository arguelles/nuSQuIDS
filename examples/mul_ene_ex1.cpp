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
#include "nuSQUIDS.h"

/*
 * This file demonstrates how use nuSQUIDS to propage neutrinos
 * considering oscillation as well as noncoherent interactions
 * such as CC/NC interactions and tau regeneration.
 */

using namespace nusquids;

int main()
{
  squids::Const units;
  const unsigned int numneu = 4;
  nuSQUIDS nus(1.*units.GeV,1.e2*units.GeV,200,numneu,neutrino,true,false);

  double phi = acos(-1.);
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track_atm = std::make_shared<EarthAtm::Track>(phi);

  nus.Set_Body(earth_atm);
  nus.Set_Track(track_atm);

  // set mixing angles and masses
  nus.Set_MixingAngle(0,1,0.563942);
  nus.Set_MixingAngle(0,2,0.154085);
  nus.Set_MixingAngle(1,2,0.785398);

  nus.Set_SquareMassDifference(1,7.65e-05);
  nus.Set_SquareMassDifference(2,0.00247);

  // sterile neutrino parameters
  nus.Set_SquareMassDifference(3,3.0);
  nus.Set_MixingAngle(1,3,0.7);

  // setup integration settings
  nus.Set_h_max( 500.0*units.km );
  //nus.Set_h_min( 1.0*nus.units.meter );
  nus.Set_rel_error(1.0e-9);
  nus.Set_abs_error(1.0e-9);

  nus.Set_ProgressBar(true);

  // construct the initial state
  marray<double,1> E_range = nus.GetERange();

  marray<double,2> inistate{E_range.size(),numneu};
  double N0 = 1.0e18;
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        // initialze muon state
        inistate[i][k] = (k == 1) ? N0*pow(E_range[i],-2) : 0.0;
      }
  }
  // set the initial state
  nus.Set_initial_state(inistate,flavor);

  nus.EvolveState();
  nus.WriteStateHDF5("./mulene_example1.hdf5");

  return 0;
}
