#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace nusquids;

int main(){
  squids::Const units;

  nuSQUIDS nus(logspace(1.e2*units.GeV,1.e6*units.GeV,59),3,both,true);

  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track_atm = std::make_shared<EarthAtm::Track>(acos(-1.));

  nus.Set_Body(earth_atm);
  nus.Set_Track(track_atm);

  marray<double,3> inistate{60,2,3};
  for(size_t i = 0 ; i < inistate.extent(0); i++){
    for(size_t r = 0; r < inistate.extent(1); r++){
      for(size_t k = 0; k < inistate.extent(2); k++){
        // initialze muon state
        inistate[i][r][k] = (k==1) ? 1:0;
      }
    }
  }
  // set the initial state
  nus.Set_initial_state(inistate,flavor);

  nus.EvolveState();

  track_atm->ReverseTrack();

  nus.EvolveState();

  for(size_t ie=0; ie<nus.GetNumE(); ie++){
    for(size_t ir = 0; ir < nus.GetNumRho(); ir++){
      for(size_t iflv=0; iflv<nus.GetNumNeu(); iflv++){
        double eps = fabs(inistate[ie][ir][iflv] - nus.EvalFlavorAtNode(iflv,ie,ir));
        if(eps>5.0e-3) // half a percent error
          std::cout << ie << " " <<  ir << " " << iflv << " " << eps << inistate[ie][ir][iflv] << " " << std::endl;
      }
    }
  }

  return 0;
}
