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

#include <memory>
#include <vector>
#include <iostream>
#include <fstream>
#include <nuSQuIDS/nuSQuIDS.h>
#include "SolarModel.h"
#include "SolarProbabilities.h"

using namespace nusquids;

int main()
{
  squids::Const units;

  std::shared_ptr<SolarModel> sm = std::make_shared<SolarModel>("Standard");
  SOP solar;
  solar.SetSolarModel(sm);

  auto e_range = linspace(0.1*units.MeV,15.*units.MeV,100);
  //auto e_range = linspace(0.1,15.,100);

  for(double enu : e_range){
    std::cout << enu/units.MeV << " " << solar.SolarOscillationProbability(enu,0.1) << " " << solar.PeeSquare(enu) << std::endl;
  }

  return 0;
}
