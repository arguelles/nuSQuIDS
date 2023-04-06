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
#include <fstream>
#include <nuSQuIDS/nuSQuIDS.h>
#include "NSI.h"

/*
 * This file demonstrates how use nuSQUIDS to propagate neutrinos
 * considering oscillation as well as non-coherent interactions
 * such as CC/NC interactions and tau regeneration for the Multiple energy case.
 */

using namespace nusquids;

int main()
{

  double big = -0.48;
  double small = -0.56;
  double DMP=-4.1666e-7;

  squids::Const units;
  unsigned int numneu=7;
  nuSQUIDSNSI nus(logspace(1.0*pow(10,small)*units.GeV,1.0*pow(10,big)*units.GeV,200),numneu,neutrino,false);
  //Here we define the trajectory that the particle follows and the object for more examples

  const double layer_1 = .49*units.km;
  std::shared_ptr<Vacuum> constdens_env0 = std::make_shared<Vacuum>();
  std::shared_ptr<Vacuum::Track> track_env0 = std::make_shared<Vacuum::Track>(layer_1);
  std::shared_ptr<ConstantDensity> constdens_env1 = std::make_shared<ConstantDensity>(1.5, 0.5); // density [gr/cm^3], ye [dimensionless]
  std::shared_ptr<ConstantDensity::Track> track_env1 = std::make_shared<ConstantDensity::Track>(layer_1);

  nus.Set_Body(constdens_env1);
  nus.Set_Track(track_env1);


    // IN order to mix the vanilla sterile and quasi-sterile neutrinos correctly we order the neutrinos as QS1, QS2, QS3, electron, muon, tau, vanilla
  // mass differences are written in the list with respect to the electron neutrino mass in the above order
  double DMDif[numneu] = {241.672, 250.005 + 7.65e-5, 258.338 + 0.00247, 1, 0,
 7.65e-5, 0.00247};

  // mixing angles are written as QS1-E, QS2-M, QS3-T, E-M, E-T, M-T, VS-E
  double DMAng[numneu] = {0.00305128,0.003,0.00295122,0.3,0.563942,0.154085,0.785398};

  for(int i=1; i<numneu; i++){
  nus.Set_SquareMassDifference(i,DMDif[i]-DMDif[0]);
  }
  //Set Mixing Angles in the way that will best emulate the paper's maths

  nus.Set_MixingAngle(0,4,DMAng[0]);
  nus.Set_MixingAngle(1,5,DMAng[1]);
  nus.Set_MixingAngle(2,6,DMAng[2]);
  nus.Set_MixingAngle(3,4,DMAng[3]);
  nus.Set_MixingAngle(4,5,DMAng[4]);
  nus.Set_MixingAngle(4,6,DMAng[5]);
  nus.Set_MixingAngle(5,6,DMAng[6]);


  std::cout << "End: setting mixing angles." << std::endl;
 
  
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
  //Se a power low spectra for the muon neutrinos (k==5), other flavors to 0.0
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        inistate[i][k] = (k == 5) ? N0*pow(E_range[i],-2) : 0.0;
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
  std::ofstream file("DMpmns.txt");
  //number of energies we want the result, notice that this can be larger than the number of the internal grid of 
  //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
  //and vacuum oscillations are solved analytically for the given energy.
  int Nen =100000;
  double lEmin=small;
  double lEmax=big;

  file << "# log10(E) E flux_i fluxRatio_i . . . ." << std::endl;
  for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
    double E=pow(10.0,lE)*units.GeV;
    file << lE << " " << E << " ";
    for(int fl=0; fl<numneu; fl++){
      file << " " <<  nus.EvalFlavor(fl, E) << " " <<  nus.EvalFlavor(fl, E)/(N0*pow(E,-2));
    }
    file << std::endl;
  }
  file.close();
  std::string plt;
  std::cout << std::endl <<  "Done! " << std::endl <<
    "  Do you want to run the gnuplot script? yes/no" << std::endl;
  std::cin >> plt;
  if(plt=="yes" || plt=="y")
    return system("./plot.plt");

  return 0;
}
