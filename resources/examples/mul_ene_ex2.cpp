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
  // We must first create a nuSQUIDS object. In order to do this
  // we must specify :
  // 1) Emin : minimum energy in GeV
  // 2) Emax : maximum energy in GeV
  // 3) Number of energy sites (integer)
  // 4) Number of neutrino flavors (integer)
  // 5) If we are going to consider neutrinos or antineutrinos or both.
  // 6) Energy log scale : true (log) or false (linear).
  // 7) Interactions : true or false.
  // We do this in the following declaration.
  nuSQUIDS nus(1.e2,1.e6,60,3,"neutrino",true,true);

  // now we need to specify a medium where this neutrinos propagate.
  // a classical scenario where this is relevant is high energy
  // atmospheric neutrinos. so we such the earth in the atmospheric
  // configuration. and we such the totally through going trayectory.

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

  nus.Set(DELTA1,0.0);

  // setup integration settings
  nus.Set_h_max( 500.0*nus.units.km );
  nus.Set_rel_error(1.0e-12);
  nus.Set_abs_error(1.0e-12);

  std::vector<double> E_range = nus.GetERange();

  // construct the initial state
  array2D inistate(60);
  double N0 = 1.0e18;
  for ( int i = 0 ; i < inistate.size(); i++){
      inistate[i].resize(3);
      for ( int k = 0; k < 3; k ++){
        // initialze muon state
        inistate[i][k] = (k == 1) ? N0*pow(E_range[i],-1.0) : 0.0;
      }
  }

  // set the initial state
  nus.Set_initial_state(inistate,"flavor");

  // we can save the current state in HDF5 format
  // for future use.
  std::cout << "Propagating nu_mu flux." << std::endl;

  nus.WriteStateHDF5("./in_muon.hdf5");
  nus.EvolveState();
  nus.WriteStateHDF5("./out_muon.hdf5");

  return 0;
}
