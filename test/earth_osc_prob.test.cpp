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

#include <iostream>
#include <iomanip>
#include <nuSQuIDS/nuSQuIDS.h>

using namespace nusquids;


void exercise_se_mode(unsigned int numneu,NeutrinoType NT, std::shared_ptr<Body> body, std::shared_ptr<Track> track){
  nuSQUIDS nus(numneu,NT);
  nus.Set_Track(track);
  nus.Set_Body(body);

  nus.Set_rel_error(1.0e-14);
  nus.Set_abs_error(1.0e-14);
  nus.Set_Basis(interaction);
  nus.Set_h(nus.units.km);
  nus.Set_h_max(300.0*nus.units.km);

  switch (numneu){
    case 3:
      nus.Set_MixingAngle(0,1,0.583996);
      nus.Set_MixingAngle(0,2,0.148190);
      nus.Set_MixingAngle(1,2,0.737324);
      nus.Set_SquareMassDifference(1,7.5e-05);
      nus.Set_SquareMassDifference(2,0.00257);
      nus.Set_CPPhase(0,2,1.);
      break;
    case 4:
      // random values for non standart parameters
      nus.Set_MixingAngle(0,1,0.583996);
      nus.Set_MixingAngle(0,2,0.148190);
      nus.Set_MixingAngle(1,2,0.737324);
      nus.Set_MixingAngle(0,3,0.1245);
      nus.Set_MixingAngle(1,3,0.5454);
      nus.Set_MixingAngle(2,3,0.32974);
      nus.Set_SquareMassDifference(1,7.5e-05);
      nus.Set_SquareMassDifference(2,0.00257);
      nus.Set_SquareMassDifference(3,1.9234);
      nus.Set_CPPhase(0,2,1.);
      nus.Set_CPPhase(0,3,0.135);
      break;
  }

  // energies
  std::vector<double> test_energies {1.0e-2,1.0e-1,1.0e0,1.0e1,1.0e2,1.0e3,1.0e4};

  std::cout << std::setprecision(3);
  std::cout << std::scientific;
  for(int flv = 0; flv < numneu; flv++){
    marray<double,1> ini_state{numneu};
    for (int iflv = 0; iflv < numneu; iflv++)
      ini_state[iflv] = ( iflv==flv ? 1.0 : 0.0 );
      for(double Enu : test_energies){
        nus.Set_E(Enu*nus.units.GeV);
        nus.Set_initial_state(ini_state,flavor);
        nus.EvolveState();
        std::cout << body->GetName() << " " << flv << " [flv] " << Enu << " [GeV] ";
        for (int i = 0; i < numneu; i++){
          double p = nus.EvalFlavor(i);
          if ( p < 1.0e-8)
            std::cout << 0.0 << " ";
          else
            std::cout << p << " ";
        }
        std::cout << std::endl;
      }
  }
}

int main(){
  // this test checks the neutrino oscillation probability
  // for 3 and 4 neutrino flavors for neutrinos/antineutrinos/both
  // and the single and multiple energy modes

  // units
  squids::Const units;

  // earth
  std::shared_ptr<Earth> earth = std::make_shared<Earth>();
  std::shared_ptr<Earth::Track> earth_track = std::make_shared<Earth::Track>(1000.0*units.km);
  //std::shared_ptr<Earth::Track> earth_track = std::make_shared<Earth::Track>(4000.0*units.km);

  exercise_se_mode(3,neutrino,earth,earth_track);
  exercise_se_mode(3,antineutrino,earth,earth_track);
  exercise_se_mode(4,neutrino,earth,earth_track);
  exercise_se_mode(4,antineutrino,earth,earth_track);

  return 0;
}
