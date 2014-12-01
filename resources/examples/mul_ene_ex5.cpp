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


class nuSQUIDSNSI: public nuSQUIDS {
  private:
    SU_vector NSI;
    std::vector<SU_vector> NSI_evol;
  public:

    void AddToPreDerive(double x){
      for(int ei = 0; ei < ne; ei++){
        SU_vector h0 = H0(E_range[ei]);
        NSI_evol[ei] = NSI.SUEvolve(h0,tunit*(x-t_ini));
      }
    }

    SU_vector HI(int ei){
      double ye = body->ye(track);
      double density = body->density(track);

      double CC = params.sqrt2*params.GF*params.Na*pow(params.cm,-3)*density*ye;
      double NC;

      if (ye < 1.0e-10){
        NC = params.sqrt2*params.GF*params.Na*pow(params.cm,-3)*density;
      }
      else {
        NC = CC*(-0.5*(1.0-ye)/ye);
      }

      // construct potential in flavor basis
      SU_vector potential = (CC+NC)*evol_b1_proj[index_rho][0][ei];
      potential += (NC)*(evol_b1_proj[index_rho][1][ei]);
      potential += (NC)*(evol_b1_proj[index_rho][2][ei]);
      // plus sign so that the NSI potential has the same sign as the VCC potential
      // and the factor of 3 comes from average n_n/n_e at Earth.
      potential += (3.0*CC)*NSI_evol[ei];

      if ((index_rho == 0 and NT=="both") or NT=="neutrino"){
          // neutrino potential
          return potential;
      } else if ((index_rho == 1 and NT=="both") or NT=="antineutrino"){
          // antineutrino potential
          return (-1.0)*potential;
      } else{
          throw std::runtime_error("nuSQUIDS::HI : unknown particle or antiparticle");
      }
    }

    nuSQUIDSNSI(double Emin,double Emax,int Esize,int numneu,std::string NT,
         bool elogscale,bool iinteraction) : nuSQUIDS(Emin,Emax,Esize,numneu,NT,elogscale,iinteraction)
    {
       assert(numneu == 3);
       // defining a complex matrix M which will contain our flavor
       // violating flavor structure.
       gsl_matrix_complex * M = gsl_matrix_complex_calloc(3,3);
       double epsilon_mutau = 1.0e-2;
       gsl_complex c { epsilon_mutau , 0.0 };
       gsl_matrix_complex_set(M,2,1,c);
       gsl_matrix_complex_set(M,1,2,gsl_complex_conjugate(c));

       NSI = SU_vector(M);

       Set(TH12,0.563942);
       Set(TH13,0.154085);
       Set(TH23,0.785398);

       // rotate to mass reprentation
       NSI.RotateToB1(params);
       NSI_evol.resize(ne);
       for(int ei = 0; ei < ne; ei++){
         NSI_evol[ei] = SU_vector(nsun);
       }
       gsl_matrix_complex_free(M);
    }
};

int main()
{
  nuSQUIDSNSI nus(1.e1,1.e3,200,3,"antineutrino",true,false);

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
  nus.Set_h_max( 200.0*nus.units.km );
  nus.Set_rel_error(1.0e-15);
  nus.Set_abs_error(1.0e-15);

  std::vector<double> E_range = nus.GetERange();

  // construct the initial state
  array2D inistate(200);
  double N0 = 1.0;
  for ( int i = 0 ; i < inistate.size(); i++){
      inistate[i].resize(3);
      for ( int k = 0; k < 3; k ++){
        // initialze muon state
        inistate[i][k] = (k == 1) ? N0 : 0.0;
      }
  }

  // set the initial state
  nus.Set_initial_state(inistate,"flavor");

  nus.Set_ProgressBar(true);
  nus.EvolveState();
  // we can save the current state in HDF5 format
  // for future use.
  nus.WriteStateHDF5("./mul_ene_ex5.hdf5");

  return 0;
}
