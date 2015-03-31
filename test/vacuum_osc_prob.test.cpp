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
#include <nuSQuIDS/nuSQUIDS.h>

using namespace nusquids;


void exercise_se_mode(unsigned int numneu,NeutrinoType NT){
  nuSQUIDS nus(numneu,NT);

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

  std::shared_ptr<Vacuum> vacuum = std::make_shared<Vacuum>();
  nus.Set_Body(vacuum);

  std::vector<double> test_baseline {1.0e-4,1.0e-3,1.0e-2,1.0e-1,1.0e0,1.0e1,1.0e2,1.0e3,1.0e4};
  std::vector<double> test_energies {1.0e-6,1.0e-5,1.0e-4,1.0e-3,1.0e-2,1.0e-1,1.0e0,1.0e1,1.0e2,1.0e3,1.0e4,1.0e5,1.0e6};

  std::cout << std::setprecision(3);
  std::cout << std::scientific;
  for(int flv = 0; flv < numneu; flv++){
    marray<double,1> ini_state{numneu};
    for (int iflv = 0; iflv < numneu; iflv++)
      ini_state[iflv] = ( iflv==flv ? 1.0 : 0.0 );
    for(double baseline :test_baseline){
      std::shared_ptr<Vacuum::Track> track_vac = std::make_shared<Vacuum::Track>(0.0,baseline*nus.units.km);
      nus.Set_Track(track_vac);
      for(double Enu : test_energies){
        nus.Set_E(Enu*nus.units.GeV);
        nus.Set_initial_state(ini_state,flavor);
        nus.EvolveState();
        std::cout << flv << " [flv] " << Enu << " [GeV] " << baseline << " [km] ";
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
}

int main(){
  // this test checks the neutrino oscillation probability
  // for 3 and 4 neutrino flavors for neutrinos/antineutrinos/both
  // and the single and multiple energy modes
  exercise_se_mode(3,neutrino);
  exercise_se_mode(3,antineutrino);
  exercise_se_mode(4,neutrino);
  exercise_se_mode(4,antineutrino);
}
