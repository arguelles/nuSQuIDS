#include <iostream>
#include <iomanip>
#include <nuSQuIDS/nuSQUIDS.h>
#include <nuSQuIDS/tools.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using namespace nusquids;

squids::Const g_units;
double density = 5.0; // gr/cm^3
double ye = 0.3;
double baseline = 700;

double ConstantDensityOscillationFormulae(NeutrinoType NT, unsigned int a, unsigned int b,
                                          double E, double L, double density_, double ye_,
                                          const gsl_matrix_complex * U, const gsl_vector * dm2){
  unsigned int numneu = U->size1;

  gsl_matrix_complex * M2 = gsl_matrix_complex_calloc(numneu,numneu);
  for(unsigned int k=1;k<numneu;k++){
    gsl_matrix_complex_set(M2,k,k,gsl_complex_rect(gsl_vector_get(dm2,k-1),0.));
  }

  gsl_matrix_complex * H = gsl_matrix_complex_alloc(numneu,numneu);
  gsl_matrix_complex * T1 = gsl_matrix_complex_alloc(numneu,numneu);

  gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,gsl_complex_rect(1.0,0.0),M2,U,gsl_complex_rect(0.0,0.0),T1);
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),U,T1,gsl_complex_rect(0.0,0.0),H);
  gsl_matrix_complex_scale(H,gsl_complex_rect(1./(2.*E),0));

  gsl_matrix_complex * V = gsl_matrix_complex_calloc(numneu,numneu);
  double CC = g_units.sqrt2*g_units.GF*g_units.Na*pow(g_units.cm,-3.)*density_*ye_;
  double NC = CC*(-0.5*(1.0-ye_)/ye_);
  for(unsigned int k=0;k<3;k++){
    if (k == 0)
      gsl_matrix_complex_set(V,k,k,gsl_complex_rect(CC+NC,0.));
    else
      gsl_matrix_complex_set(V,k,k,gsl_complex_rect(NC,0.));
  }
  if (NT==antineutrino)
    gsl_matrix_complex_scale(V,gsl_complex_rect(-1.,0.));

  gsl_matrix_complex_add(H,V);

  gsl_matrix_complex_scale(H,gsl_complex_rect(0.,-L));
  //std::cout << "H" << std::endl;
  //gsl_matrix_complex_print(H);

  gsl_matrix_complex * eH = gsl_matrix_complex_alloc(numneu,numneu);
  squids::math_detail::matrix_exponential(eH,H);
  //std::cout << "eH" << std::endl;
  //gsl_matrix_complex_print(eH);

  gsl_vector_complex * inistate = gsl_vector_complex_calloc(numneu);
  gsl_vector_complex_set(inistate,a,GSL_COMPLEX_ONE);
  gsl_vector_complex * outstate = gsl_vector_complex_calloc(numneu);
  gsl_blas_zgemv(CblasNoTrans,GSL_COMPLEX_ONE,eH,inistate,GSL_COMPLEX_ZERO,outstate);

  double p = gsl_complex_abs2(gsl_vector_complex_get(outstate,b));

  gsl_vector_complex_free(inistate);
  gsl_vector_complex_free(outstate);
  gsl_matrix_complex_free(V);
  gsl_matrix_complex_free(eH);
  gsl_matrix_complex_free(T1);
  gsl_matrix_complex_free(M2);
  gsl_matrix_complex_free(H);

  return p;
}

void exercise_se_mode(unsigned int numneu,NeutrinoType NT, std::shared_ptr<Body> body, std::shared_ptr<Track> track){
  nuSQUIDS nus(numneu,NT);
  squids::Const units;

  nus.Set_Track(track);
  nus.Set_Body(body);

  nus.Set_rel_error(1.0e-15);
  nus.Set_abs_error(1.0e-15);
  nus.Set_h(units.km);
  nus.Set_h_max(300.0*units.km);

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
  //std::vector<double> test_energies {1.0e-2,1.0e-1,1.0e0,1.0e1,1.0e2,1.0e3,1.0e4};
  std::vector<double> test_energies {1.0e0};

  // set up simple formulae parameters
  auto U = nus.GetTransformationMatrix();
  if(NT==antineutrino){
    gsl_matrix_complex_conjugate(U.get());
  }
  gsl_vector * dm2 = gsl_vector_alloc(numneu-1);
  for(unsigned int ii = 1; ii < numneu; ii++){
    gsl_vector_set(dm2,ii-1,nus.Get_SquareMassDifference(ii));
  }

  std::cout << std::setprecision(3);
  std::cout << std::scientific;
  for(int iflv = 0; iflv < numneu; iflv++){
    marray<double,1> ini_state{numneu};
    for (int flv = 0; flv < numneu; flv++)
      ini_state[flv] = ( iflv==flv ? 1.0 : 0.0 );

    for(double Enu : test_energies){
      nus.Set_E(Enu*units.GeV);
      nus.Set_initial_state(ini_state,flavor);
      nus.EvolveState();
      for (unsigned int fflv = 0; fflv < numneu; fflv++){
          double p_nsq = nus.EvalFlavor(fflv);
          double p_for = ConstantDensityOscillationFormulae(NT,iflv,fflv,Enu*units.GeV,baseline*units.km,density,ye,U.get(),dm2);
          if ( std::abs(p_nsq - p_for) > 1.0e-4 )
            std::cout << NT << " " <<  iflv << " " << fflv << " " << Enu << " " << baseline << " "<< p_nsq << " " << p_for << std::endl;
        }
    }
  }
}

int main(){
  // this test checks the neutrino oscillation probability
  // for 3 and 4 neutrino flavors for neutrinos/antineutrinos/both
  // and the single and multiple energy modes

  // units
  squids::Const units;

  // constant density
  std::shared_ptr<ConstantDensity> constdens = std::make_shared<ConstantDensity>(density,ye);
  std::shared_ptr<ConstantDensity::Track> track_constdens = std::make_shared<ConstantDensity::Track>(baseline*units.km);

  exercise_se_mode(3,neutrino,constdens,track_constdens);
  exercise_se_mode(3,antineutrino,constdens,track_constdens);
  exercise_se_mode(4,neutrino,constdens,track_constdens);
  exercise_se_mode(4,antineutrino,constdens,track_constdens);

  return 0;
}
