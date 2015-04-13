#include <nuSQuIDS/nuSQUIDS.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace nusquids;

int main(){

  nuSQUIDS nus1(1.e2,1.e6,60,3,neutrino,true,false);
  auto e_range = nus1.GetERange();

  nuSQUIDS nus2(e_range,3,neutrino,false);
  for (unsigned int ie = 0; ie < nus1.GetNumE(); ie++){
    double delta = std::abs(nus1.GetERange()[ie] -  nus2.GetERange()[ie]);
    if (delta > 1.0e-15)
      std::cout << ie << " " << delta << std::endl;
  }
  // i need to put some things inside in order to write it out
  double phi = acos(-0.5);
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track_atm = std::make_shared<EarthAtm::Track>(phi);
  nus1.Set_Body(earth_atm);
  nus1.Set_Track(track_atm);
  marray<double,2> inistate{60,3};
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        // initialze muon state
        inistate[i][k] = 1.;
      }
  }
  // set the initial state
  nus1.Set_initial_state(inistate,flavor);

  nus1.WriteStateHDF5("mult_energy_test.hdf5");
  nuSQUIDS nus3("mult_energy_test.hdf5");
  for (unsigned int ie = 0; ie < nus1.GetNumE(); ie++){
    double delta = std::abs(nus1.GetERange()[ie] -  nus3.GetERange()[ie]);
    if (delta > 1.0e-15)
      std::cout << ie << " " << delta << std::endl;
  }

  nuSQUIDS nus4(std::move(nus2));
  for (unsigned int ie = 0; ie < nus1.GetNumE(); ie++){
    double delta = std::abs(nus1.GetERange()[ie] -  nus4.GetERange()[ie]);
    if (delta > 1.0e-15)
      std::cout << ie << " " << delta << std::endl;
  }

  return 0;
}
