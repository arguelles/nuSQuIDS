#include "SolarProbabilities.h"

using namespace squids;
using namespace nusquids;

double SOP::SolarOscillationProbability(double E,double r) const {
  //if (params == nullptr)
  //  throw std::runtime_error("No oscillation parameters set");
  if (solar_model == nullptr)
    throw std::runtime_error("No solar model set");

  SU_vector h = Hamiltonian(E,r);
  h.RotateToB0(params);
  auto eigensyst = h.GetEigenSystem();
  // order according to eigenvalues
  gsl_eigen_hermv_sort(eigensyst.first.get(),eigensyst.second.get(),GSL_EIGEN_SORT_VAL_ASC);

  auto UPMNS = GetTransformationMatrix();

  double osc_prob = 0;
  for(unsigned int i = 0; i < numneu; i++){
    osc_prob += gsl_complex_abs2(gsl_matrix_complex_get(eigensyst.second.get(),nue,i))*\
                gsl_complex_abs2(gsl_matrix_complex_get(UPMNS.get(),nue,i));
  }
  return osc_prob;
}

SU_vector SOP::Hamiltonian(double E, double r) const {

  double electron_number_density = solar_model->eDensity(r);
  SU_vector H = DM2*(1./(2.*E));
  double CC_prefactor = HI_constants;// defined in nusquids
  double CC = CC_prefactor*electron_number_density;
  H += CC*b1_proj[neutrino][nue];

  return std::move(H);
}

double SOP::PeeSquare(double E) const{
  double normalization = 1./RadialIntegratedFluxes(E);
  return integrate([&](double r){
              double sum_flux = 0.;
              for ( unsigned int i = 0; i < solar_model->NumComp(); i++){
                if (not solar_model->isline[SolarModel::FluxType(i)] ) {
                  sum_flux += SolarOscillationProbability(E,r)*solar_model->nuFlux(r,E,SolarModel::FluxType(i));
                } else if ( std::abs(E - solar_model->spectrum_limits[SolarModel::FluxType(i)][0])/E < 1.0e-2 ) {
                  // if close to the line evaluate the line
                  double E_line = solar_model->spectrum_limits[SolarModel::FluxType(i)][0];
                  sum_flux += SolarOscillationProbability(E_line,r)*solar_model->nuFlux(r,E_line,SolarModel::FluxType(i));
                }
              }
              return sum_flux;
          },0.,0.99,1e-3)*normalization;
}

double SOP::RadialIntegratedFluxes(double E) const{
  return integrate([&](double r){
              double sum_flux = 0.;
              for ( unsigned int i = 0; i < solar_model->NumComp(); i++){
                if (not solar_model->isline[SolarModel::FluxType(i)] ) {
                  sum_flux += solar_model->nuFlux(r,E,SolarModel::FluxType(i));
                } else if ( std::abs(E - solar_model->spectrum_limits[SolarModel::FluxType(i)][0])/E < 1.0e-2 ) {
                  // if close to the line evaluate the line
                  sum_flux += solar_model->nuFlux(r,solar_model->spectrum_limits[SolarModel::FluxType(i)][0],SolarModel::FluxType(i));
                }
              }
              return sum_flux;
          },0.,0.99,1e-3);
}
