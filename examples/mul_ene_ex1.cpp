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
  nuSQUIDS nus(1.,1.e2,200,4,nuSQUIDS::neutrino,true,false);

  double phi = acos(-1.);
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track_atm = std::make_shared<EarthAtm::Track>(phi);

  nus.Set_Body(earth_atm);
  nus.Set_Track(track_atm);

  // set mixing angles and masses
  nus.Set(TH12,0.563942);
  nus.Set(TH13,0.154085);
  nus.Set(TH23,0.785398);

  nus.Set(DM21SQ,7.65e-05);
  nus.Set(DM31SQ,0.00247);

  // sterile neutrino parameters
  nus.Set(DM41SQ, 3.0);
  nus.Set(TH24, 0.7);

  nus.Set(DELTA1,0.0);

  // setup integration settings
  nus.Set_h_max( 500.0*nus.units.km );
  //nus.Set_h_min( 1.0*nus.units.meter );
  nus.Set_rel_error(1.0e-9);
  nus.Set_abs_error(1.0e-9);

  // construct the initial state
  marray<double,1> E_range = nus.GetERange();

  marray<double,2> inistate{E_range.size(),3};
  double N0 = 1.0e18;
  for ( int i = 0 ; i < inistate.size(); i++){
      for ( int k = 0; k < 3; k ++){
        // initialze muon state
        inistate[i][k] = (k == 1) ? N0*pow(E_range[i],-2) : 0.0;
      }
  }
  // set the initial state
  nus.Set_initial_state(inistate,"flavor");

  nus.EvolveState();
  nus.WriteStateHDF5("./mulene_example1.hdf5");

  return 0;
}
