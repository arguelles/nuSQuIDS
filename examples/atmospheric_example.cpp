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
 * This file demostrate the use of nuSQUIDSAtm class
 * which is design to handle a bundle of trayectories
 * such as those involved in atmospheric neutrino
 * experiments.
 */

using namespace nusquids;

int main()
{
  int numneu = 3;
  nuSQUIDSAtm nus_atm(-1.,0.2,20,5.e2,1.e6,150,numneu,nuSQUIDS::both,true,true);

  // set mixing angles and masses
  nus_atm.Set(TH12,0.563942);
  nus_atm.Set(TH13,0.154085);
  nus_atm.Set(TH23,0.785398);

  nus_atm.Set(DM21SQ,7.65e-05);
  nus_atm.Set(DM31SQ,0.00247);

  nus_atm.Set(DELTA1,0.0);

  if( numneu > 3){
    nus_atm.Set(DM41SQ,1.0);
    nus_atm.Set(TH24,0.26237);
  }

  // setup integration settings
  nus_atm.Set_rel_error(1.0e-11);
  nus_atm.Set_abs_error(1.0e-11);

  // construct the initial state
  double N0 = 1.0;
  array4D inistate(nus_atm.GetNumCos());
  for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
    inistate[ci].resize(nus_atm.GetNumE());
    for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
      inistate[ci][ei].resize(2);
      for ( int rho = 0; rho < 2; rho ++ ){
      inistate[ci][ei][rho].resize(3);
        for (int flv = 0; flv < numneu; flv++){
          // initialze muon state
          inistate[ci][ei][rho][flv] = (flv == 1) ? N0 : 0.0;
        }
      }
    }
  }

  // set the initial state
  nus_atm.Set_initial_state(inistate,"flavor");

  nus_atm.Set_ProgressBar(true);
  nus_atm.EvolveState();
  // we can save the current state in HDF5 format
  // for future use.
  nus_atm.WriteStateHDF5("./atmospheric_example_numneu_"+std::to_string(numneu)+".hdf5");

  return 0;
}
