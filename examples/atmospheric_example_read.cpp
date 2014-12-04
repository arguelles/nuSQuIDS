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
#include <fstream>
#include "nuSQUIDS.h"

/*
 * This file demostrate the use of nuSQUIDSAtm class
 * which is design to handle a bundle of trayectories
 * such as those involved in atmospheric neutrino
 * experiments.
 */

using namespace nusquids;

int main(int argc,  char** argv)
{
  if(argc != 2){
    std::cerr << "Argument number not valid! Given: " <<  argc << std::endl;
    exit(1);
  }

  std::string filename = (std::string) argv[1];

  nuSQUIDSAtm nus_atm(filename);

  double costh = -0.9;
  int nu_e = 0;
  int nu_mu = 1;
  int nu_tau = 2;
  int neutrino = 0;

  std::ofstream output("outputfile.dat");

  for(double log10E = 3.0; log10E < 6.0; log10E +=  0.1 ){
    double enu = pow(10.0,log10E);
    output << enu << " ";
    output << nus_atm.EvalFlavor(nu_e,costh,enu,neutrino) << " ";
    output << nus_atm.EvalFlavor(nu_mu,costh,enu,neutrino) << " ";
    output << nus_atm.EvalFlavor(nu_tau,costh,enu,neutrino) << " ";
    output << std::endl;
  }

  return 0;
}
