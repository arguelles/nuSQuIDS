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


class nuSQUIDSLV: public nuSQUIDS {
  //friend class nuSQUIDS;
  public:
    SU_vector HI(int ei){
      //std::cout << 1.0e-21*E_range[ei]*evol_b1_proj[index_rho][1][ei] << std::endl;
      return 1.0e-21*E_range[ei]*evol_b1_proj[index_rho][1][ei];
    }

    nuSQUIDSLV(double Emin_,double Emax_,int Esize_,int numneu_,string NT_,
         bool elogscale_,bool iinteraction_)
    {
       Init(Emin_,Emax_,Esize_,numneu_,NT_,
           elogscale_,iinteraction_);
    }
};

int main()
{
  nuSQUIDSLV nus(1.e2,1.e6,150,3,"neutrino",true,false);

  double phi = acos(-1.);
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track_atm = std::make_shared<EarthAtm::Track>(phi);

  nus.Set_Body(earth_atm);
  nus.Set_Track(track_atm);

  // set mixing angles and masses
  nus.Set("th12",0.563942);
  nus.Set("th13",0.154085);
  nus.Set("th23",0.785398);

  nus.Set("dm21sq",7.65e-05);
  nus.Set("dm31sq",0.00247);

  nus.Set("delta1",0.0);

  // setup integration settings
  nus.Set("h_max", 500.0*nus.units.km );
  nus.Set("rel_error", 1.0e-15);
  nus.Set("abs_error", 1.0e-15);

  vector<double> E_range = nus.GetERange();

  // construct the initial state
  array2D inistate(150);
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

  //nus.WriteStateHDF5("./mul_ene_ex4_in.hdf5");
  nus.EvolveState();
  nus.WriteStateHDF5("./mul_ene_ex4.hdf5");

  return 0;
}
