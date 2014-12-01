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
 * This file demostrates how to read a nuSQUIDS file
 * and evalute its contents.
 */

using namespace nusquids;

int main()
{
  nuSQUIDS nus("/Users/carguelles/Workspace/SQuIDS/git_version/nuSQuIDS/mul_ene_ex5.hdf5");
  nus.WriteStateHDF5("/Users/carguelles/Workspace/SQuIDS/git_version/nuSQuIDS/test.hdf5");

  nus.WriteStateHDF5("/Users/carguelles/Workspace/SQuIDS/git_version/nuSQuIDS/test.hdf5", "track1");
  //nuSQUIDS nus2("/Users/carguelles/Workspace/SQuIDS/git_version/nuSQuIDS/test.hdf5","track1");

  return 0;
}

