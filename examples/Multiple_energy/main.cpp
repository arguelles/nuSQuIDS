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
 * This file demonstrates how use nuSQUIDS to propagate neutrinos
 * considering oscillation as well as non-coherent interactions
 * such as CC/NC interactions and tau regeneration for the Multiple energy case.
 */

using namespace nusquids;

int main()
{
  squids::Const units;
  //Number of neutrinos, 4 is the case with one sterile neutrino
  std::cout << "(3) Three Active Neutrinos, " << "(4) 3+1 Three Active and One Sterile Neutrino" << std::endl;
  unsigned int numneu;
  std::cin >>numneu;
  if( not(numneu==3 || numneu==4)){
    throw std::runtime_error("Only (3) or (4) are valid options");
  }
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

  if(numneu==4){
    // sterile neutrino parameters
    nus.Set_SquareMassDifference(3,0.1);
    nus.Set_MixingAngle(1,3,0.1);
  }
  
  //Here we set the maximum size for the integration step, important for fast or sharp variations of the density.
  nus.Set_h_max( 500.0*units.km );

  //Setting the numerical precision of gsl integrator.
  nus.Set_rel_error(1.0e-10);
  nus.Set_abs_error(1.0e-10);

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
  //Propagate the neutrinos in the earth for the path defined in path
  nus.EvolveState();
  //This functions save the stat of the system, in this case after the evolution
  std::cout << std::endl << "writing the outputs..." << std::endl;


  //In this part we will save the values in a txt file to be able to plot or manipulate later.
  //Notice that this is not going to have all the information about the quantum evolution, for that 
  //we need to save the information using the HDF5 writing function.
  std::ofstream file("fluxes_flavor.txt");
  //number of energies we want the result, notice that this can be larger than the number of the internal grid of 
  //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
  //and vacuum oscillations are solved analytically for the given energy.
  int Nen =1000;
  double lEmin=0;
  double lEmax=4;
  
  file << "# log10(E) E flux_i fluxRatio_i . . . ." << std::endl;
  for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
    double E=pow(10.0,lE)*units.GeV;
    file << lE << " " << E << " ";
    for(int fl=0; fl<numneu; fl++){
      file << " " <<  nus.EvalFlavor(fl, E) << " " <<  nus.EvalFlavor(fl, E)/(N0*pow(E,-2));
    }
    file << std::endl;
  }
  std::string plt;
  std::cout << std::endl <<  "Done! " << std::endl <<
    "  Do you want to run the gnuplot script? yes/no" << std::endl;
  std::cin >> plt;
  if(plt=="yes" || plt=="y")
    return system("./plot.plt");

  return 0;
}
