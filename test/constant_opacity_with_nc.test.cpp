#include <iostream>
#include <iomanip>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/tools.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_odeiv2.h>

using namespace nusquids;
squids::Const units;

const double density = 5.0; // gr/cm^3
const double ye = 0.3;
const double baseline = 12000;
// here the factor of two is because we are averaging over neutrons and protons mass. The number of nucleons
// should be the sum of neutrons and protons.
const double num_nucleon = density*(units.gr*pow(units.cm,-3))*(2.0/(units.proton_mass+units.neutron_mass));

double flux_function(double enu){
  return 1.;
}

struct GSLParamsHelper {
  double * delE;
  size_t number_of_energy_nodes;
  nuSQUIDS::InteractionStructure * is;
};

int DiffEquationKernel(double t,
                       const double * f, double * dfdE,
                       void *params){
  (void)(t); /* avoid unused parameter warning */
  GSLParamsHelper* cast_params = (GSLParamsHelper*)(params);
  double * delE = cast_params->delE;
  unsigned int number_of_energy_nodes = cast_params->number_of_energy_nodes;
  nuSQUIDS::InteractionStructure* is = cast_params->is;
  for(unsigned int ei = 0; ei < number_of_energy_nodes; ei++)
    dfdE[ei] = -f[ei]*num_nucleon*(is->sigma_CC[0][0][ei]+is->sigma_NC[0][0][ei]);

  for(unsigned int e2 = 1; e2 < number_of_energy_nodes; e2++){
    for(unsigned int e1 = 0; e1 < e2; e1++){
      // this only does neutrinos
      dfdE[e1] += f[e2]*delE[e2-1]*is->dNdE_NC[0][0][e2][e1]*is->sigma_NC[0][0][e2]*num_nucleon;
    }
  }

  return GSL_SUCCESS;
}

std::vector<double> SimplePropagateFlux(marray<double,1> energy_nodes,
                                        std::shared_ptr<nuSQUIDS::InteractionStructure> intstruct,
                                        double tf){
  std::vector<double> phi(energy_nodes.size());
  std::vector<double> delE(energy_nodes.size()-1);
  for(unsigned int ei=0; ei < phi.size(); ei++){
    phi[ei] =  flux_function(energy_nodes[ei]);
  }
  for(unsigned int ei=0; ei < delE.size(); ei++)
    delE[ei] = energy_nodes[ei+1]-energy_nodes[ei];

  GSLParamsHelper parameters {delE.data(),energy_nodes.size(),intstruct.get()};
  gsl_odeiv2_system sys = {DiffEquationKernel,nullptr,energy_nodes.size(),&parameters};

  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,1.0e-12,1.0e-12,0.0);
  double t=0;
  int gsl_status = gsl_odeiv2_driver_apply(d, &t,tf, phi.data());
  gsl_odeiv2_driver_free(d);

  if( gsl_status != GSL_SUCCESS ){
    throw std::runtime_error("ERROR:TEST GSL Evolution failed.");
  }

  return phi;
}

int main(){
  // units
  const unsigned int numneu = 3;
  double Emin = 1.0e3*units.GeV;
  double Emax = 1.0e6*units.GeV;
  unsigned int Esize = 100;
  nuSQUIDS nus(Emin,Emax,Esize,numneu,both,true,true);

  // constant density
  std::shared_ptr<ConstantDensity> constdens = std::make_shared<ConstantDensity>(density,ye);
  std::shared_ptr<ConstantDensity::Track> track_constdens = std::make_shared<ConstantDensity::Track>(baseline*units.km);

  // setting the track and body
  nus.Set_Body(constdens);
  nus.Set_Track(track_constdens);

  auto e_range = nus.GetERange();
  marray<double,3> inistate{nus.GetNumE(),2,numneu};
  std::fill(inistate.begin(),inistate.end(),0);
  for ( unsigned int ei = 0 ; ei < nus.GetNumE(); ei++){
    for ( unsigned int rho = 0; rho < 2; rho ++ ){
      for ( unsigned int flv = 0; flv < numneu; flv++){
        // initialze muon state
        inistate[ei][rho][flv] = flux_function(e_range[ei]);
      }
    }
  }

  nus.Set_rel_error(1.0e-20);
  nus.Set_abs_error(1.0e-20);

  nus.Set_initial_state(inistate,flavor);
  nus.Set_IncludeOscillations(false);
  //nus.Set_IncludeOscillations(true);
  nus.Set_OtherRhoTerms(true);
  nus.EvolveState();

  auto int_structure = nus.GetInteractionStructure();

  std::vector<double> phi = SimplePropagateFlux(e_range,int_structure,track_constdens->GetFinalX());

  for(unsigned int ei = 0; ei < nus.GetNumE(); ei++){
    double phi_e = nus.EvalFlavorAtNode(0,ei,0);
    double phi_mu = nus.EvalFlavorAtNode(1,ei,0);
    double phi_tau = nus.EvalFlavorAtNode(2,ei,0);

    if( std::abs(phi_e - phi[ei])/phi[ei] > 5.0e-2 or
        std::abs(phi_mu - phi[ei])/phi[ei] > 5.0e-2 or
        std::abs(phi_tau - phi[ei])/phi[ei] > 5.0e-2){
      if(!(std::abs(phi_e - phi[ei])<1.0e-2))
        std::cout << e_range[ei] << " " <<  phi_e << " " << phi_mu << " " << phi_tau << " " << phi[ei] << std::endl;
    }
  }

  return 0;
}
