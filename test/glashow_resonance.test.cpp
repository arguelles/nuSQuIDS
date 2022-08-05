#include <nuSQuIDS/xsections.h>
#include <nuSQuIDS/nuSQuIDS.h>
#include <SQuIDS/const.h>
#include <SQuIDS/SUNalg.h>

std::ostream& operator<<(std::ostream &f, const std::vector<double> &vec)
{
  f << "[";
  for (const double &v : vec)
    f << " " << v;
  f << "]";
  return f;
}

double integrate_xs(double start, double stop, double step, double e0)
{
  using namespace nusquids;
  if(start>stop){
    std::cout << "Integration stop is less than start." << std::endl;
    return 0;
  }

  unsigned int nsteps=(stop-start)/step;
  double xs = 0;
  const squids::Const constants;
  GlashowResonanceCrossSection gr;
  for (unsigned int i=0; i < nsteps; i++) {
    double logE = start + step*i;
    double dE = pow(10., logE+step) - pow(10., logE);
    xs += gr.SingleDifferentialCrossSection(e0,pow(10., logE)*constants.GeV,
                                            NeutrinoCrossSections::electron,
                                            NeutrinoCrossSections::antineutrino,
                                            NeutrinoCrossSections::GR)*dE;
  }
  double scale = xs/gr.TotalCrossSection(e0,
   NeutrinoCrossSections::electron,
   NeutrinoCrossSections::antineutrino,
   NeutrinoCrossSections::GR)/gr.WDecayBranchingFraction(NeutrinoCrossSections::muon);
  return scale;
}

int main (int argc, char const *argv[])
{
  using namespace nusquids;

  const squids::Const constants;
  std::cout << "Checking single-differential cross-section" << std::endl;
  std::cout << "\\int dsigma/dE / sigma = " << integrate_xs(0, 7, 0.001, 6.3e6*constants.GeV) << std::endl;
  assert(std::abs(integrate_xs(0, 7, 0.001, 6.3e6*constants.GeV) - 1) < 3e-3);

  std::cout << "Evolving a nue_bar line spectrum at the Glashow resonance" << std::endl;
  const unsigned int numneu = 3;
  const unsigned int num_steps = 201;
  squids::Const units;
  std::shared_ptr<NullCrossSections> ndcs = std::make_shared<NullCrossSections>();
  nuSQUIDS squid(logspace(1e4*units.GeV,1e7*units.GeV,num_steps),numneu,both,true,ndcs);

  std::shared_ptr<ConstantDensity> const_dens = std::make_shared<ConstantDensity>(10.,0.5);
  std::shared_ptr<ConstantDensity::Track> const_dens_track = std::make_shared<ConstantDensity::Track>(10000.*units.km);

  squid.Set_Body(const_dens);
  squid.Set_Track(const_dens_track);

  squid.Set_IncludeOscillations(false);
  squid.Set_GlashowResonance(true);
  squid.Set_PositivityConstrain(false);
  squid.Set_TauRegeneration(false);

  // setup integration settings
  squid.Set_h_max( 100.0*units.km );
  squid.Set_GSL_step(gsl_odeiv2_step_rkf45);

  squid.Set_rel_error(1.0e-25);
  squid.Set_abs_error(1.0e-25);

  // construct the initial state
  marray<double,3> inistate{num_steps,2,numneu};
  std::fill(inistate.begin(), inistate.end(), 0.);
  marray<double,1> E_range = squid.GetERange();
  for(unsigned int ie = 0; ie < E_range.size(); ie++){
    inistate[ie][1][0] = 1./(E_range[ie]);
  }

  // set the initial state
  squid.Set_initial_state(inistate,flavor);

  std::cout << "Propagating. . . " << std::endl;
  squid.EvolveState();
  std::cout << "Propagation done" << std::endl;
  std::ostream &output = std::cout;
  for (auto i=0; i < num_steps-1; i++) {
    double ee = E_range[i];
    double de = (E_range[i+1]-E_range[i])/units.GeV;
    output << i << " ";
    output << E_range[i]/units.GeV << " ";
    output << squid.EvalFlavorAtNode(0,i,0)*ee << " ";
    output << squid.EvalFlavorAtNode(1,i,0)*ee << " ";
    output << squid.EvalFlavorAtNode(2,i,0)*ee << " ";
    output << squid.EvalFlavorAtNode(0,i,1)*ee << " ";
    output << squid.EvalFlavorAtNode(1,i,1)*ee << " ";
    output << squid.EvalFlavorAtNode(2,i,1)*ee << " ";
    output << std::endl;
    /*
    for (auto flav=0; flav<3; flav++) {
      assert(squid.EvalFlavorAtNode(flav,i,0) == 0);
      // Extra interaction channel suppresses nue_bar more than other flavors
      //assert(squid.EvalFlavorAtNode(flav,i,1) >= squid.EvalFlavorAtNode(0,i,1));
      if(squid.EvalFlavorAtNode(flav,i,1) < squid.EvalFlavorAtNode(0,i,1))
        output << "nu_e bar not smaller than other flavors: " << flav << ' ' << i << ' '
         << squid.EvalFlavorAtNode(flav,i,1) << ' ' << squid.EvalFlavorAtNode(0,i,1) << std::endl;;
      // If there is some nue_bar flux left over, there must also be some other flavors
      if (flav!=0 && squid.EvalFlavorAtNode(0,i,1)>0 && squid.EvalFlavorAtNode(flav,i,1)<=0){
        output << "Unexpected lack of flux: " << flav << ' ' << i
         << ' ' << squid.EvalFlavorAtNode(0,i,1) << ' ' << squid.EvalFlavorAtNode(flav,i,1) << std::endl;
      }
    }
    */
  }

  return 0;
}
