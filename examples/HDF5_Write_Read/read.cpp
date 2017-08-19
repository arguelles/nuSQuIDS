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

/*
 * This file demonstrates how use nuSQUIDS can read the full state of the system
 * from an hdf5 file and use the information to safe a test file with the final fluxes.
 */

using namespace nusquids;

int main()
{
  squids::Const units;

  //Here we create the nusquids object reading the state from the hdf5 file.
  //nuSQUIDS inus("./initial_state.hdf5");
  //nuSQUIDS fnus("./final_state.hdf5");

  nuSQUIDSAtm<> inus("./pion_atmospheric_initial.hdf5");
  nuSQUIDSAtm<> fnus("./pion_atmospheric_2441_1.000000_0.000000_0.160875_0.000000_0.000000_0.000000.hdf5");
  
  //In this part we will save the values in a txt file to be able to plot or manipulate later.
  //Notice that this is not going to have all the information about the quantum evolution, for that 
  //we need to save the information using the HDF5 writing function.
  std::ofstream file("fluxes_flavor.txt");
  //number of energies we want the result, notice that this can be larger than the number of the internal grid of 
  //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
  //and vacuum oscillations are solved analytically for the given energy.
  unsigned int Nen =1000;
  double lEmin=0;
  double lEmax=4;
  unsigned int numneu=inus.GetNumNeu();  

  file << "# log10(E) E flux_i fluxRatio_i . . . ." << std::endl;
  for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
    double E=pow(10.0,lE)*units.GeV;
    file << lE << " " << E << " ";
    for(int fl=0; fl<numneu; fl++){
      //the ration in the second column is with the initial muon flux, the others are zero.
      file << " " <<  fnus.EvalFlavor(fl, -1,E) << " " <<  fnus.EvalFlavor(fl,-1, E)/inus.EvalFlavor(1, -1,E);
    }
    file << std::endl;
  }
  std::string plt;
  std::cout << std::endl <<  "Done! " << std::endl <<
    "  Do you want to run the gnuplot script? yes/no" << std::endl;
  std::cin >> plt;
  if(plt=="yes" || plt=="y")
    return system("./plot.plt");

}
