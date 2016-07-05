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
 * This file demonstrates how use nuSQUIDS can write the full state of the system, 
 * is the same computation as in the multiple energy mode but now storing the full
 * initial and final state in a file.
 */

using namespace nusquids;

int main()
{
  squids::Const units;

  //Number of neutrinos
  unsigned int numneu=3;

  //Declaration of the nuSQUIDS object, the arguments are:
  //(1). Minimum energy
  //(2). Maximum energy
  //(3). Number of energy bins
  //(4). Number of neutrino states
  //(5). Neutrino or anti-neutrino case
  //(6). Energy logarithmic scale? 
  //(7). Scattering non coherent interactions. 
  nuSQUIDS nus(1.*units.GeV,1.e4*units.GeV,200,numneu,neutrino,true,false);
  
  //Here we define the trajectory that the particle follows and the object for more examples
  // of how construct a track and object look body_track example.
  //zenith angle, neutrinos crossing the earth
  double phi = acos(-1.);
  //Declaration of the body, EarthAtm is one of the predefined bodies
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  //Definition of the track, in encodes the trajectory inside the body, here is declared with the zenith angle.
  std::shared_ptr<EarthAtm::Track> track_atm = std::make_shared<EarthAtm::Track>(phi);
  //We set this in the nusSQuID object.
  nus.Set_Body(earth_atm);
  nus.Set_Track(track_atm);


  // set mixing angles and masses
  nus.Set_MixingAngle(0,1,0.563942);
  nus.Set_MixingAngle(0,2,0.154085);
  nus.Set_MixingAngle(1,2,0.785398);
  nus.Set_SquareMassDifference(1,7.65e-05);
  nus.Set_SquareMassDifference(2,0.00247);

#ifdef STERILE
  // sterile neutrino parameters
  nus.Set_SquareMassDifference(3,0.1);
  nus.Set_MixingAngle(1,3,0.1);
#endif

  //Here we set the maximum size for the integration step, important for fast or sharp variations of the density.
  nus.Set_h_max( 500.0*units.km );

  //We set the GSL step function
  nus.Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus.Set_rel_error(1.0e-5);
  nus.Set_abs_error(1.0e-5);

  //Set true the progress bar during the evolution.
  nus.Set_ProgressBar(true);

  //Construct the initial state
  //E_range is an array that contains all the energies.
  marray<double,1> E_range = nus.GetERange();
  //Array that contains the initial state of the system, fist component is energy and second every one of the flavors
  marray<double,2> inistate{E_range.size(),numneu};
  double N0 = 1.0e18;
  //Se a power low spectra for the muon neutrinos (k==1), other flavors to 0.0
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        inistate[i][k] = (k == 1) ? N0*pow(E_range[i],-2) : 0.0;
      }
  }

  //Set the initial state in nuSQuIDS object
  nus.Set_initial_state(inistate,flavor);
  nus.WriteStateHDF5("./initial_state.hdf5");
  //Propagate the neutrinos in the earth for the path defined in path
  nus.EvolveState();

  //This functions save the stat of the system, in this case after the evolution
  nus.WriteStateHDF5("./final_state.hdf5");

  return 0;
}
