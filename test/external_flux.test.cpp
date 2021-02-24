#include <iostream>
#include <iomanip>
#include <vector>

#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/resources.h>

using namespace nusquids;

// Exponential source
// The assumption will be that it will be a all flavor neutrino source
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
  if(interaction_length>=std::numeric_limits<double>::max()){
    //when lambda_1 -> infinity the phi(t) term can be ignored
    return decay_length*(1-exp(-t/decay_length));
  }
  double c=(interaction_length*decay_length)/(decay_length-interaction_length);
  return -c*exp(-t/interaction_length)+c*exp(-t/decay_length);
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
      if(std::abs(rel) > 1.e-5 or std::isnan(del))
        std::cout << fl << ' ' << del << '\n';
    }
  }
}

int main(){
  squids::Const units;
  const double integration_tol=1e-8;

  {// interactions true and osc true, vac
    std::cout << "Interactions true, oscillations true, vacuum" << std::endl;
    nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,neutrino,true);
    double decay_length = 100.0*units.meter;
    std::shared_ptr<EmittingSlab> emit = std::make_shared<EmittingSlab>(decay_length,0.0,0.5);
    std::shared_ptr<EmittingSlab::Track> track = std::make_shared<EmittingSlab::Track>(1.0*units.km);
    nus.Set_Body(emit);
    nus.Set_Track(track);
    nus.Set_rel_error(integration_tol);
    nus.Set_abs_error(integration_tol);
    nus.Set_NeutrinoSources(true);
    nus.Set_IncludeOscillations(true);

    marray<double,1> E_range = nus.GetERange();
    marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
    std::fill(inistate.begin(),inistate.end(),0.0);
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    test(nus,1000,2,6,std::numeric_limits<double>::max(),decay_length,track);
  }

  {// interactions true and osc false, vac
    std::cout << "Interactions true, oscillations false, vacuum" << std::endl;
    nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,neutrino,true);
    double decay_length = 100.0*units.meter;
    std::shared_ptr<EmittingSlab> emit = std::make_shared<EmittingSlab>(decay_length,0.0,0.5);
    std::shared_ptr<EmittingSlab::Track> track = std::make_shared<EmittingSlab::Track>(1.0*units.km);
    nus.Set_Body(emit);
    nus.Set_Track(track);
    nus.Set_rel_error(integration_tol);
    nus.Set_abs_error(integration_tol);
    nus.Set_NeutrinoSources(true);
    nus.Set_IncludeOscillations(false);

    marray<double,1> E_range = nus.GetERange();
    marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
    std::fill(inistate.begin(),inistate.end(),0.0);
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    test(nus,1000,2,6,std::numeric_limits<double>::max(),decay_length,track);
  }

  {// interactions false  and osc false, vac
    std::cout << "Interactions false, oscillations false, vacuum" << std::endl;
    nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,neutrino,false);
    double decay_length = 100.0*units.meter;
    std::shared_ptr<EmittingSlab> emit = std::make_shared<EmittingSlab>(decay_length,0.0,0.5);
    std::shared_ptr<EmittingSlab::Track> track = std::make_shared<EmittingSlab::Track>(1.0*units.km);
    nus.Set_Body(emit);
    nus.Set_Track(track);
    nus.Set_rel_error(integration_tol);
    nus.Set_abs_error(integration_tol);
    nus.Set_NeutrinoSources(true);
    nus.Set_IncludeOscillations(false);

    marray<double,1> E_range = nus.GetERange();
    marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
    std::fill(inistate.begin(),inistate.end(),0.0);
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    test(nus,1000,2,6,std::numeric_limits<double>::max(),decay_length,track);
  }

  {// interactions false  and osc true, vac
    std::cout << "Interactions false, oscillations true, vacuum" << std::endl;
    nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,neutrino,false);
    double decay_length = 100.0*units.meter;
    std::shared_ptr<EmittingSlab> emit = std::make_shared<EmittingSlab>(decay_length,0.0,0.5);
    std::shared_ptr<EmittingSlab::Track> track = std::make_shared<EmittingSlab::Track>(1.0*units.km);
    nus.Set_Body(emit);
    nus.Set_Track(track);
    nus.Set_rel_error(integration_tol);
    nus.Set_abs_error(integration_tol);
    nus.Set_NeutrinoSources(true);
    nus.Set_IncludeOscillations(true);

    marray<double,1> E_range = nus.GetERange();
    marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
    std::fill(inistate.begin(),inistate.end(),0.0);
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    test(nus,1000,2,6,std::numeric_limits<double>::max(),decay_length,track);
  }

  // We can't write down a convenient analytical solution to the differential 
  // equation for the neutrino flux when NC cascading is taking effect, so to
  // test that we do the right thing in the presence of absorption we use a 
  // trick cross section which is only non-zero for charged current interactions
  // (which purely remove neutrinos in our treatment)
  class CCOnlyCrossSections : public nusquids::NeutrinoDISCrossSectionsFromTables{
    using Parent=nusquids::NeutrinoDISCrossSectionsFromTables;
  public:
    CCOnlyCrossSections(std::string p):Parent(p){}
    double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override{
      if(current!=CC)
        return 0;
      return Parent::TotalCrossSection(Enu, flavor, neutype, current);
    }
    double SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override{
      if(current!=CC)
        return 0;
      return Parent::SingleDifferentialCrossSection(E1, E2, flavor, neutype, current);
    }
  };
  std::shared_ptr<CrossSectionLibrary> ccOnlyXS(new CrossSectionLibrary);
  {
    std::string xsdir = nusquids::getResourcePath()+"/xsections/";
    ccOnlyXS->addTarget(nusquids::proton, CCOnlyCrossSections(xsdir+"csms_proton.h5"));
    ccOnlyXS->addTarget(nusquids::neutron,CCOnlyCrossSections(xsdir+"csms_neutron.h5"));
  }
  
  {// interactions true and osc false, matter
    std::cout << "Interactions true, oscillations false, matter" << std::endl;
    nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,neutrino,true);
    nus.SetNeutrinoCrossSections(ccOnlyXS);
    double decay_length = 100.0*units.meter;
    double density = 5.0;
    //use a material which has no neutrons to reduce the number of cross sections we need to deal with
    std::shared_ptr<EmittingSlab> emit = std::make_shared<EmittingSlab>(decay_length,density,1.0);
    std::shared_ptr<EmittingSlab::Track> track = std::make_shared<EmittingSlab::Track>(100*units.km);
    nus.Set_Body(emit);
    nus.Set_Track(track);
    nus.Set_rel_error(integration_tol);
    nus.Set_abs_error(integration_tol);
    nus.Set_NeutrinoSources(true);
    nus.Set_IncludeOscillations(false);

    marray<double,1> E_range = nus.GetERange();
    marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
    std::fill(inistate.begin(),inistate.end(),0.0);
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    auto int_struct = nus.GetInteractionStructure();
    double target_number_density = density*units.gr*pow(units.cm,-3)/(units.proton_mass + units.electron_mass);
    
    for(unsigned int ie = 0; ie < E_range.size(); ie++){
      double E=E_range[ie];
      for(unsigned int fl=0; fl<nus.GetNumNeu(); fl++){
        //we used only proton targets and CC cross sections
        double invlen = int_struct->sigma_CC[0][0][fl][ie]*target_number_density;
        double interaction_length = 1./invlen;
        double exact = analytical_solution(track->GetFinalX(),interaction_length,decay_length);
        double computed = nus.EvalFlavor(fl, E); 
        double del = computed - exact;
        double rel = del/exact;
        if(std::abs(rel) > 1.e-5 or std::isnan(del))
          std::cout << ie << ' ' << fl << ' ' << del <<  ' ' << rel << ' ' << exact << ' ' << computed << '\n';
      }
    }
  }
  
  {// interactions true and osc true, matter
    std::cout << "Interactions true, oscillations true, matter" << std::endl;
    nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,neutrino,true);
    nus.SetNeutrinoCrossSections(ccOnlyXS);
    double decay_length = 100.0*units.meter;
    double density = 5.0;
    //use a material which has no neutrons to reduce the number of cross sections we need to deal with
    std::shared_ptr<EmittingSlab> emit = std::make_shared<EmittingSlab>(decay_length,density,1.0);
    std::shared_ptr<EmittingSlab::Track> track = std::make_shared<EmittingSlab::Track>(100*units.km);
    nus.Set_Body(emit);
    nus.Set_Track(track);
    nus.Set_rel_error(integration_tol);
    nus.Set_abs_error(integration_tol);
    nus.Set_NeutrinoSources(true);
    nus.Set_IncludeOscillations(true);

    marray<double,1> E_range = nus.GetERange();
    marray<double,2> inistate{E_range.size(),nus.GetNumNeu()};
    std::fill(inistate.begin(),inistate.end(),0.0);
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    auto int_struct = nus.GetInteractionStructure();
    double target_number_density = density*units.gr*pow(units.cm,-3)/(units.proton_mass + units.electron_mass);
    
    for(unsigned int ie = 0; ie < E_range.size(); ie++){
      double E=E_range[ie];
      for(unsigned int fl=0; fl<nus.GetNumNeu(); fl++){
        //we used only proton targets and CC cross sections
        double invlen = int_struct->sigma_CC[0][0][fl][ie]*target_number_density;
        double interaction_length = 1./invlen;
        double exact = analytical_solution(track->GetFinalX(),interaction_length,decay_length);
        double computed = nus.EvalFlavor(fl, E); 
        double del = computed - exact;
        double rel = del/exact;
        if(std::abs(rel) > 1.e-5 or std::isnan(del))
          std::cout << ie << ' ' << fl << ' ' << del <<  ' ' << rel << ' ' << exact << ' ' << computed << '\n';
      }
    }
  }

  return 0;
}
