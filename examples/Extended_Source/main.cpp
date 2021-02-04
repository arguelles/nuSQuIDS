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
#include "nu_source.h"

/*
 * This file demonstrates how use nuSQUIDS to propagate neutrinos
 * considering oscillation as well as non-coherent interactions
 * such as CC/NC interactions and tau regeneration for the Multiple energy case.
 */

using namespace nusquids;

int main()
{
  squids::Const units;
  nuSQUIDS nus(logspace(1.*units.GeV,1.e4*units.GeV,200),3,neutrino,true);
  // Defined the parameters of our emittor
  double decay_length = 100.0*units.meter;
  unsigned int flv = 1; // nu mu
  //Body that we have defined in the header
  std::shared_ptr<EmittingVacuum> emit_vac = std::make_shared<EmittingVacuum>(flv,decay_length);
  //Since our body inherits from Vacuum, we can use Vacuum::Track instead of defining an specific EmittingVacuum::Track
  //This works because our inheritted body has the same geometry as the parent body.
  std::shared_ptr<EmittingVacuum::Track> track_vac = std::make_shared<EmittingVacuum::Track>(1.0*units.km);
  //We set this in the nuSQuIDS object.
  nus.Set_Body(emit_vac);
  nus.Set_Track(track_vac);

  // set mixing angles and masses
  nus.Set_MixingAngle(0,1,0.563942);
  nus.Set_MixingAngle(0,2,0.154085);
  nus.Set_MixingAngle(1,2,0.785398);
  nus.Set_SquareMassDifference(1,7.65e-05);
  nus.Set_SquareMassDifference(2,0.00247);

  //We set the GSL step function
  nus.Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus.Set_rel_error(1.0e-5);
  nus.Set_abs_error(1.0e-5);

  //Set true the progress bar during the evolution.
  nus.Set_ProgressBar(true);

  //Enable having neutrinos emitted from thbe body.
  nus.Set_NeutrinoSources(true);

  //Construct the initial state
  //E_range is an array that contains all the energies.
  marray<double,1> E_range = nus.GetERange();
  //Array that contains the initial state of the system, fist component is energy and second every one of the flavors
  marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
  std::fill(inistate.begin(),inistate.end(),0.0);

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
    for(int fl=0; fl<nus.GetNumNeu(); fl++){
      file << " " <<  nus.EvalFlavor(fl, E) << " " <<  nus.EvalFlavor(fl, E);
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
