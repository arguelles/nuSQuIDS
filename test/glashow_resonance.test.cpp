

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
  
  int nsteps=(stop-start)/step;
  double xs = 0;
  const squids::Const constants;
  GlashowResonanceCrossSection gr;
  for (int i=0; i < nsteps; i++) {
    double logE = start + step*i;
    double dE = pow(10., logE+step) - pow(10., logE);
    xs += gr.SingleDifferentialCrossSection(e0,pow(10., logE)*constants.GeV,
NeutrinoCrossSections::electron,NeutrinoCrossSections::antineutrino,NeutrinoCrossSections::GR)*dE;
  }
  double scale = xs/gr.TotalCrossSection(e0,
      NeutrinoCrossSections::electron,
      NeutrinoCrossSections::antineutrino,
      NeutrinoCrossSections::GR)/gr.GetMuonicBranchingFraction();
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
  const unsigned int num_steps = 151;
  squids::Const units;
  nuSQUIDS squid(1e4*units.GeV,1e7*units.GeV,num_steps,numneu,both,true,true);

  double phi = acos(-1);
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track_atm = std::make_shared<EarthAtm::Track>(phi);
  
  squid.Set_Body(earth_atm);
  squid.Set_Track(track_atm);
  
  // set mixing angles and masses
  squid.Set_MixingAngle(0,1,0.563942);
  squid.Set_MixingAngle(0,2,0.154085);
  squid.Set_MixingAngle(1,2,0.785398);

  squid.Set_SquareMassDifference(1,7.65e-05);
  squid.Set_SquareMassDifference(2,0.00247);

  // setup integration settings
  squid.Set_h_max( 500.0*units.km );
  squid.Set_rel_error(1.0e-24);
  squid.Set_abs_error(1.0e-24);
  
  // construct the initial state
  marray<double,3> inistate{num_steps,2,numneu};
  std::fill(inistate.begin(), inistate.end(), 0.);
  marray<double,1> E_range = squid.GetERange();
  int e0 = 140;
  inistate[e0][1][0] = 1./(E_range[e0+1]-E_range[e0])*units.GeV;

  // set the initial state
  squid.Set_initial_state(inistate,flavor);
  squid.Set_IncludeOscillations(false);
  squid.Set_PositivityConstrain(false);
  
  squid.EvolveState();
  std::ostream &output = std::cout;
  for (auto i=0; i < num_steps-1; i++) {
    double de = (E_range[i+1]-E_range[i])/units.GeV;
    output << i << " ";
    output << E_range[i]/units.GeV << " ";
    output << squid.EvalFlavorAtNode(0,i,0)*de << " ";
    output << squid.EvalFlavorAtNode(1,i,0)*de << " ";
    output << squid.EvalFlavorAtNode(2,i,0)*de << " ";
    output << squid.EvalFlavorAtNode(0,i,1)*de << " ";
    output << squid.EvalFlavorAtNode(1,i,1)*de << " ";
    output << squid.EvalFlavorAtNode(2,i,1)*de << " ";
    output << std::endl;
    for (auto flav=0; flav<3; flav++) {
      assert(squid.EvalFlavorAtNode(flav,i,0) == 0);
      // Extra interaction channel suppresses nue_bar more than other flavors
      assert(squid.EvalFlavorAtNode(flav,i,1) >= squid.EvalFlavorAtNode(0,i,1));
      // If there is some nue_bar flux left over, there must also be some other flavors
      if (squid.EvalFlavorAtNode(0,i,1) != 0)
        assert(squid.EvalFlavorAtNode(flav,i,1) > 0);
    }
  }

  return 0;
}