#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace nusquids;

// Exponential source
// The assumption will be that it will be a single flavor neutrino source
// with an exponential decay profile.
class EmittingSlab: public ConstantDensity{
  private:
    const double decay_length;
  public :
    EmittingSlab(double decay_length, double density, double ye):
      ConstantDensity(density,ye),decay_length(decay_length){}
    void injected_neutrino_flux(marray<double,3>& flux, const GenericTrack& track, const nuSQUIDS& nusquids) override {
      double x_cur = track.GetX();
      for(unsigned int ei=0; ei < nusquids.GetNumE(); ei++){
        for(unsigned int rhoi = 0; rhoi < nusquids.GetNumRho(); rhoi++){
          for(unsigned int flv = 0; flv < nusquids.GetNumNeu(); flv++){
            //flux[ei][rhoi][flv] = (flv == flavor) ? exp(-x_cur/decay_length) : 0.0;
            flux[ei][rhoi][flv] = exp(-x_cur/decay_length);
          }
        }
      }
    }
};

// Analytical solution
// Solution of a differential equation of the form
// dphi(t)/dt = - phi(t)/lambda_1 + exp(-t/lambda_2)
// with phi(0) = 0 and where
// lambda_1 : interaction_length
// lambda_2 : decay_length
double analytical_solution(double t, double interaction_length, double decay_length){
  double c = 1./decay_length - 1./interaction_length;
  return (1 - exp(-t*c))/c;
}

void test(const nuSQUIDS& nus,
           unsigned int Nen, double lEmin, double lEmax,double interaction_length, double decay_length,
           std::shared_ptr<EmittingSlab::Track> track){
  squids::Const units;
  for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
    double E=pow(10.0,lE)*units.GeV;
    for(int fl=0; fl<nus.GetNumNeu(); fl++){
      double exact = analytical_solution(track->GetFinalX(),interaction_length,decay_length);
      double del = nus.EvalFlavor(fl, E) - exact;
      double rel = del/exact;
      //std::cout << fl << " " << exact << " " << nus.EvalFlavor(fl, E) << " " << del << std::endl;
      if(fabs(rel) > 1.e-5 or std::isnan(del))
        std::cout << fl << " " << del << std::endl;
    }
  }
}

int main(){
  squids::Const units;

  {// interactions true and osc true, vac
    nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,neutrino,true);
    double decay_length = 100.0*units.meter;
    //unsigned int flv = 1; // nu mu
    std::shared_ptr<EmittingSlab> emit = std::make_shared<EmittingSlab>(decay_length,0.0,0.5);
    std::shared_ptr<EmittingSlab::Track> track = std::make_shared<EmittingSlab::Track>(1.0*units.km);
    nus.Set_Body(emit);
    nus.Set_Track(track);
    nus.Set_rel_error(1.0e-8);
    nus.Set_abs_error(1.0e-8);
    nus.Set_GSL_step(gsl_odeiv2_step_rk4);
    nus.Set_NeutrinoSources(true);
    nus.Set_IncludeOscillations(true);
    //nus.Set_IncludeOscillations(false);
    //nus.Set_ProgressBar(true);

    marray<double,1> E_range = nus.GetERange();
    marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
    std::fill(inistate.begin(),inistate.end(),0.0);
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    test(nus,1000,2,6,std::numeric_limits<double>::max(),decay_length,track);
  }

  {// interactions true and osc false, vac
    nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,neutrino,true);
    double decay_length = 100.0*units.meter;
    //unsigned int flv = 1; // nu mu
    std::shared_ptr<EmittingSlab> emit = std::make_shared<EmittingSlab>(decay_length,0.0,0.5);
    std::shared_ptr<EmittingSlab::Track> track = std::make_shared<EmittingSlab::Track>(1.0*units.km);
    nus.Set_Body(emit);
    nus.Set_Track(track);
    nus.Set_rel_error(1.0e-8);
    nus.Set_abs_error(1.0e-8);
    nus.Set_GSL_step(gsl_odeiv2_step_rk4);
    nus.Set_NeutrinoSources(true);
    nus.Set_IncludeOscillations(false);
    //nus.Set_IncludeOscillations(false);
    //nus.Set_ProgressBar(true);

    marray<double,1> E_range = nus.GetERange();
    marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
    std::fill(inistate.begin(),inistate.end(),0.0);
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    test(nus,1000,2,6,std::numeric_limits<double>::max(),decay_length,track);
  }

  {// interactions false  and osc false, vac
    nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,neutrino,false);
    double decay_length = 100.0*units.meter;
    //unsigned int flv = 1; // nu mu
    std::shared_ptr<EmittingSlab> emit = std::make_shared<EmittingSlab>(decay_length,0.0,0.5);
    std::shared_ptr<EmittingSlab::Track> track = std::make_shared<EmittingSlab::Track>(1.0*units.km);
    nus.Set_Body(emit);
    nus.Set_Track(track);
    nus.Set_rel_error(1.0e-8);
    nus.Set_abs_error(1.0e-8);
    nus.Set_GSL_step(gsl_odeiv2_step_rk4);
    nus.Set_NeutrinoSources(true);
    nus.Set_IncludeOscillations(false);
    //nus.Set_IncludeOscillations(false);
    //nus.Set_ProgressBar(true);

    marray<double,1> E_range = nus.GetERange();
    marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
    std::fill(inistate.begin(),inistate.end(),0.0);
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    test(nus,1000,2,6,std::numeric_limits<double>::max(),decay_length,track);
  }

  {// interactions false  and osc true, vac
    nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,neutrino,false);
    double decay_length = 100.0*units.meter;
    //unsigned int flv = 1; // nu mu
    std::shared_ptr<EmittingSlab> emit = std::make_shared<EmittingSlab>(decay_length,0.0,0.5);
    std::shared_ptr<EmittingSlab::Track> track = std::make_shared<EmittingSlab::Track>(1.0*units.km);
    nus.Set_Body(emit);
    nus.Set_Track(track);
    nus.Set_rel_error(1.0e-8);
    nus.Set_abs_error(1.0e-8);
    nus.Set_GSL_step(gsl_odeiv2_step_rk4);
    nus.Set_NeutrinoSources(true);
    nus.Set_IncludeOscillations(true);
    //nus.Set_IncludeOscillations(false);
    //nus.Set_ProgressBar(true);

    marray<double,1> E_range = nus.GetERange();
    marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
    std::fill(inistate.begin(),inistate.end(),0.0);
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    test(nus,1000,2,6,std::numeric_limits<double>::max(),decay_length,track);
  }

  {// interactions true and osc true, matter
    nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,neutrino,true);
    double decay_length = 100.0*units.meter;
    //unsigned int flv = 1; // nu muo
    double density = 5.0;
    std::shared_ptr<EmittingSlab> emit = std::make_shared<EmittingSlab>(decay_length,density,1.0);
    std::shared_ptr<EmittingSlab::Track> track = std::make_shared<EmittingSlab::Track>(100.0*units.km);
    nus.Set_Body(emit);
    nus.Set_Track(track);
    nus.Set_rel_error(1.0e-10);
    nus.Set_abs_error(1.0e-10);
    nus.Set_GSL_step(gsl_odeiv2_step_rk4);
    nus.Set_NeutrinoSources(true);
    nus.Set_IncludeOscillations(false);
    //nus.Set_IncludeOscillations(false);
    //nus.Set_ProgressBar(true);

    marray<double,1> E_range = nus.GetERange();
    marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
    std::fill(inistate.begin(),inistate.end(),0.0);
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    // computing/getting the interaction length
    auto int_struct = nus.GetInteractionStructure();
    double target_number_density = 2.*density*units.gr*pow(units.cm,-3)/(units.proton_mass + units.neutron_mass + units.electron_mass);
    /*
    std::cout << int_struct->sigma_NC.extent(0) << std::endl;
    std::cout << int_struct->sigma_NC[0][0][0][10] << " " << int_struct->sigma_NC[1][0][0][10] << std::endl;
    std::cout << int_struct->sigma_CC[0][0][0][10] << " " << int_struct->sigma_CC[1][0][0][10] << std::endl;
    return 0;
    */
    for(unsigned int ie = 0; ie < E_range.size(); ie++){
      double E=E_range[ie];
      for(unsigned int fl=0; fl<nus.GetNumNeu(); fl++){
        double interaction_length = 1./((int_struct->sigma_NC[0][0][fl][ie] + int_struct->sigma_CC[0][0][fl][ie])*target_number_density);
        //interaction_length = std::numeric_limits<double>::max();
        double exact = analytical_solution(track->GetFinalX(),interaction_length,decay_length);
        double del = nus.EvalFlavor(fl, E) - exact;
        double rel = del/exact;
        //std::cout << fl << " " << exact << " " << nus.EvalFlavor(fl, E) << " " << del << std::endl;
        if(fabs(rel) > 1.e-5 or std::isnan(del))
          std::cout << fl << " " << del <<  " " << rel << std::endl;
      }
    }
  }

  return 0;
}
