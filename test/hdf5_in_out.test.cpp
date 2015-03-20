#include <nuSQuIDS/nuSQUIDS.h>
#include <iostream>
#include <iomanip> 
#include <vector>

using namespace nusquids;

int main(){

  std::cout << std::setprecision(3);
  std::cout << std::fixed;

  nuSQUIDS nus(1.e2,1.e6,60,3,neutrino,true,true);

  double phi = acos(-0.5);
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track_atm = std::make_shared<EarthAtm::Track>(phi);

  nus.Set_Body(earth_atm);
  nus.Set_Track(track_atm);

  // set mixing angles and masses
  nus.Set_MixingParametersToDefault();

  // setup integration settings
  nus.Set_h_max( 500.0*nus.units.km );
  nus.Set_rel_error(1.0e-12);
  nus.Set_abs_error(1.0e-12);

  marray<double,1> E_range = nus.GetERange();

  // construct the initial state
  marray<double,2> inistate{60,3};
  double N0 = 1.0e18;
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        // initialze muon state
        inistate[i][k] = N0*pow(E_range[i],-1.0);
      }
  }

  // set the initial state
  nus.Set_initial_state(inistate,flavor);

  nus.EvolveState();
  nus.WriteStateHDF5("./out_muon_test.hdf5");

  nuSQUIDS nus_read("./out_muon_test.hdf5");

  for ( int i = 0 ; i < nus_read.GetNumE(); i++){
      for ( int k = 0; k < nus_read.GetNumNeu(); k ++){
        // initialze muon state
        std::cout << fabs(nus.EvalFlavorAtNode(k,i)/inistate[i][k] - nus_read.EvalFlavorAtNode(k,i)/inistate[i][k]) << " ";
      }
      std::cout << std::endl;
  }

  return 0;
}
