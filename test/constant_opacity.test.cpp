#include <iostream>
#include <iomanip>
#include <nuSQuIDS/nuSQUIDS.h>
#include <nuSQuIDS/tools.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using namespace nusquids;

double density = 5.0; // gr/cm^3
double ye = 0.3;
double baseline = 12000;

double flux_function(double enu){
  return 1.;
}

int main(){
  // units
  const unsigned int numneu = 3;
  squids::Const units;
  double Emin = 1.0e3*units.GeV;
  double Emax = 1.0e6*units.GeV;
  unsigned int Esize = 50;
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
  nus.Set_OtherRhoTerms(false);
  nus.EvolveState();

  auto int_structure = nus.GetInteractionStructure();
  double number_density = (units.gr*pow(units.cm,-3))*density*2.0/(units.proton_mass+units.neutron_mass);
  double column_density = baseline*units.km*number_density;

  for(unsigned int ei = 0 ; ei < nus.GetNumE(); ei++){
    for(unsigned int rho = 0; rho < 2; rho++){
      for(unsigned int flv = 0; flv < numneu; flv++){
        // initialze muon state
        double exp_analytic = exp(-column_density*(int_structure->sigma_CC[rho][flv][ei]+int_structure->sigma_NC[rho][flv][ei]));
        double exp_nusquids = nus.EvalFlavorAtNode(flv,ei,rho);
        if( std::abs(exp_nusquids - exp_analytic)/exp_analytic > 5.0e-3 )
          std::cout << flv << ' ' << rho << ' ' << ei << ' ' << exp_nusquids << " " << exp_analytic <<" " << exp_nusquids/exp_analytic << std::endl;
      }
    }
  }

  return 0;
}
