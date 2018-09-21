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
 * This file demonstrates how to the astrophysical flavo ratio by means
 * of using the averaged out approximation. We do this in two ways in this
 * example: first explicitly using the PMNS matrix and the formula given
 * in the literature and second by means of nuSQuIDS fast averaging functionality.
 * For simplicity we will do this example in the single energy mode, but it
 * can be performed in the multiple energy mode too.
 */

using namespace nusquids;

int main()
{
  // We must first create a nuSQUIDS object. In order to do this
  // we must specify the number of neutrino flavors and if we are
  // going to consider neutrino or antineutrino oscillations.

  // In this example we set N_neutrino = 3 and Type = "neutrino".
  nuSQUIDS nus(3,neutrino);

  // We will now set the mixing parameters and square mass
  // differences. The angles will be given in radiant and the
  // differences in eV^2.

  // We use the standard parametrization as described in the
  // documentation.

  // mixing angles
  nus.Set_MixingAngle(0,1,0.563942);
  nus.Set_MixingAngle(0,2,0.154085);
  nus.Set_MixingAngle(1,2,0.785398);
  // square mass differences
  nus.Set_SquareMassDifference(1,7.65e-05);
  nus.Set_SquareMassDifference(2,0.00247);
  // CP phase
  nus.Set_CPPhase(0,2,0.0);

  // Define a Const object for handleing the units
  squids::Const units;

  // Now we set the neutrino energy which we are interested on.
  // Energies are always given in natural units. To handle
  // the units the nuSQUIDS object has a unit subclass
  // which contains the most used units.

  nus.Set_E(1.0*units.PeV);

  // To calculate atmospheric neutrino oscillation probabilities
  // we need to specify a different body.


  // The neutrinos will propagate in Vacuum. So we do
  std::shared_ptr<Vacuum> vacuum = std::make_shared<Vacuum>();
  std::shared_ptr<Vacuum::Track> vacuum_track = std::make_shared<Vacuum::Track>(1.e3*units.kparsec);
  nus.Set_Body(vacuum);
  nus.Set_Track(vacuum_track);

  //Here we set the initial state for the flavor, a pion produced flavor composition
  marray<double,1> ini_state({3},{1,2,0});
  nus.Set_initial_state(ini_state,flavor);

  std::cout << "In state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  // We do the calculation
  nus.EvolveState();

  // nuSQuIDS can calculate the average oscillation probability
  // where the oscillation frequencies are larger than some value
  // we call this value scale. When this mode is used its often
  // valuable to know which frequencies have been averaged out.
  // NuSQuIDS this by modifiying a boolean vector provided to the
  // evaluation function.

  double scale = 0.;
  std::vector<bool> is_avg(3);

  // No averaging can be done by setting calling the same function
  // without the scale argument or by setting the scale argument
  // as large as possible.
  scale = std::numeric_limits<double>::max();

  // Output the result
  std::cout << "Out state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i, scale, is_avg);
      std::cout << " (" << (is_avg[i] ? "avg." : "no avg.")  << ") ";
    }
    std::cout << std::endl;
  }

  // The case of full averaging can be obtained by setting the
  // scale to zero.

  scale = 0.0;

  // Output the result
  std::cout << "Out state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i, scale, is_avg);
      std::cout << " (" << (is_avg[i] ? "avg." : "no avg.")  << ") ";
    }
    std::cout << std::endl;
  }

  // In this case we can also perform the calculation
  // by of the averaged oscillation probability
  // formulae. In this case the flavor ratio is given by
  // the initial flavor ratio and the PMNS matrix elements.

  auto PMNS = nus.GetTransformationMatrix();

  // The PMNS matrix is returned as a std::unique_ptr<gsl_matrix_complex>.
  // You can obtain the raw pointer by using the .get() member function
  // that returns a raw gsl_matrix_compplex pointer. We return a std::unique_ptr
  // instead of a raw pointer since the former deallocates itself automatically
  // thus been less error prone.

  std::cout << nus.GetERange().front()/units.GeV << " ";
  for(unsigned int flv_final = 0; flv_final  < nus.GetNumNeu(); flv_final++){
    double flv_content_final = 0.0;
    for(unsigned int imass = 0; imass < nus.GetNumNeu(); imass++){
      for(unsigned int flv_ini= 0; flv_ini < nus.GetNumNeu(); flv_ini++){
        flv_content_final += gsl_complex_abs2(gsl_matrix_complex_get(PMNS.get(),flv_final,imass))*
                             gsl_complex_abs2(gsl_matrix_complex_get(PMNS.get(),flv_ini,imass))*
                             ini_state[flv_ini];
      }
    }
    std::cout << flv_content_final << "  ";
  }
  std::cout << std::endl;

  return 0;
}
