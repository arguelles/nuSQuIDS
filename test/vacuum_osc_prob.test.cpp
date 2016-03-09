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

#define SQ(x) ((x)*(x))

using namespace nusquids;

double VacuumOscillationFormulae(NeutrinoType NT, unsigned int a, unsigned int b, double E, double L, gsl_matrix_complex * U, gsl_vector * dm2){
  double p = 0; double dm2_ij;
  gsl_complex U_product;

  if(a == b)
    p++;

  for(unsigned int j = 0; j < U->size1; j++){
    for(unsigned int i = j+1; i < U->size1; i++){
      if (j == 0)
        dm2_ij = gsl_vector_get(dm2,i-1);
      else
        dm2_ij = gsl_vector_get(dm2,i-1) - gsl_vector_get(dm2,j-1);

      if ( NT == neutrino) {
        U_product = gsl_complex_mul(
                                    gsl_complex_mul(gsl_complex_conjugate(gsl_matrix_complex_get(U,a,i)),
                                                    gsl_matrix_complex_get(U,b,i)),
                                    gsl_complex_mul(gsl_matrix_complex_get(U,a,j),
                                                    gsl_complex_conjugate(gsl_matrix_complex_get(U,b,j)))
                                   );
      } else {
        // for antineutrinos U -> U^\dagger
        U_product = gsl_complex_mul(
                                    gsl_complex_mul(gsl_complex_conjugate(gsl_matrix_complex_get(U,b,i)),
                                                    gsl_matrix_complex_get(U,a,i)),
                                    gsl_complex_mul(gsl_matrix_complex_get(U,b,j),
                                                    gsl_complex_conjugate(gsl_matrix_complex_get(U,a,j)))
                                   );
      }


      p += 2.*GSL_IMAG(U_product)*sin(dm2_ij*L/(2.*E));
      p -= 4.*GSL_REAL(U_product)*SQ(sin(dm2_ij*L/(4.*E)));
    }
  }
  return p;
}

void exercise_se_mode(unsigned int numneu, NeutrinoType NT){
  squids::Const units;
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

  // set up simple formulae parameters
  auto U = nus.GetTransformationMatrix();
  gsl_vector * dm2 = gsl_vector_alloc(numneu-1);
  for(unsigned int ii = 1; ii < numneu; ii++){
    gsl_vector_set(dm2,ii-1,nus.Get_SquareMassDifference(ii));
  }

  std::cout << std::setprecision(3);
  std::cout << std::scientific;
  marray<double,1> ini_state{numneu};

  for(unsigned int iflv = 0; iflv < numneu; iflv++){
    // setting up the initial state
    for (unsigned int ii = 0; ii < numneu; ii++)
      ini_state[ii] = ( ii==iflv ? 1.0 : 0.0 );

    // checking each baseline and energy
    for(double baseline : test_baseline){
      std::shared_ptr<Vacuum::Track> track_vac = std::make_shared<Vacuum::Track>(baseline*units.km);
      nus.Set_Track(track_vac);
      for(double Enu : test_energies){
        nus.Set_E(Enu*units.GeV);
        nus.Set_initial_state(ini_state,flavor);
        nus.EvolveState();
        for (unsigned int fflv = 0; fflv < numneu; fflv++){
          double p_nsq = nus.EvalFlavor(fflv);
          double p_for = VacuumOscillationFormulae(NT,iflv,fflv,Enu*units.GeV,baseline*units.km,U.get(),dm2);
          if ( std::abs(p_nsq - p_for) > 1.0e-4 )
            std::cout << NT << " " <<  iflv << " " << fflv << " " << Enu << " " << baseline << " "<< p_nsq << " " << p_for << std::endl;
        }
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
