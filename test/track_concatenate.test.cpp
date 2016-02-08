#include <nuSQuIDS/nuSQUIDS.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace nusquids;

int main(){

  squids::Const units;
  std::cout << std::setprecision(3);
  std::cout << std::fixed;

  nuSQUIDS nus(1.e0*units.GeV,1.e1*units.GeV,60,3,neutrino,false,false);
  const double distance = 500.*units.km;
  std::shared_ptr<Vacuum> vac = std::make_shared<Vacuum>();
  std::shared_ptr<Vacuum::Track> track = std::make_shared<Vacuum::Track>(distance);

  nus.Set_Body(vac);
  nus.Set_Track(track);

  // set mixing angles and masses
  nus.Set_MixingParametersToDefault();

  // setup integration settings
  nus.Set_h_max( 500.0*units.km );
  nus.Set_rel_error(1.0e-12);
  nus.Set_abs_error(1.0e-12);

  // construct the initial state
  marray<double,2> inistate{nus.GetNumE(),nus.GetNumNeu()};
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        inistate[i][k] = (k == 1);
      }
  }

  nus.Set_initial_state(inistate,flavor);
  nus.EvolveState();

  marray<double,2> outstate_1{nus.GetNumE(),nus.GetNumNeu()};
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        outstate_1[i][k] = nus.EvalFlavorAtNode(k,i);
      }
  }

  nus.Set_initial_state(inistate,flavor);
  track = std::make_shared<Vacuum::Track>(distance/2);
  nus.Set_Track(track);
  nus.EvolveState();
  nus.Set_Track(track);
  nus.EvolveState();

  for ( int i = 0 ; i < nus.GetNumE(); i++){
      for ( int k = 0; k < nus.GetNumNeu(); k ++){
       // std::cout << outstate_1[i][k] << " ";
       double dif = fabs(nus.EvalFlavorAtNode(k,i) - outstate_1[i][k]);
       if (dif > 1.0e-15)
        std::cout << i << " " << k << " " << dif << std::endl;
      }
  }

  return 0;
}
