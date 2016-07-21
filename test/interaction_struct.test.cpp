#include <SQuIDS/Const.h>
#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace nusquids;

int main(){
  squids::Const units;

  nuSQUIDS nus1(logspace(1.e2*units.GeV,1.e6*units.GeV,60),3,both,true);
  auto e_range = nus1.GetERange();
  nus1.InitializeInteractions();
  auto int_struct = nus1.GetInteractionStructure();

  marray<double,1> costh_range {{10},{-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1}};
  nuSQUIDSAtm<> nusatm(costh_range,e_range,3,both,int_struct);

  if(*(nusatm.GetInteractionStructure)() != (*int_struct)){
    std::cout << "Fail 1" << std::endl;
  }

  double N0 = 1.0;
  marray<double,4> inistate{nusatm.GetNumCos(),nusatm.GetNumE(),nusatm.GetNumRho(),nusatm.GetNumNeu()};
  std::fill(inistate.begin(),inistate.end(),0);

  for ( int ci = 0 ; ci < nusatm.GetNumCos(); ci++){
    for ( int ei = 0 ; ei < nusatm.GetNumE(); ei++){
      for ( int rho = 0; rho < nusatm.GetNumRho(); rho ++ ){
        for (int flv = 0; flv < nusatm.GetNumNeu(); flv++){
          // initialze muon state
          inistate[ci][ei][rho][flv] = (flv == 1) ? N0 : 0.0;
        }
      }
    }
  }

  // set the initial state
  nusatm.Set_initial_state(inistate,flavor);

  nusatm.WriteStateHDF5("atm_test.hdf5");

  nuSQUIDSAtm<> nusatm_2("atm_test.hdf5");

  if(*(nusatm_2.GetInteractionStructure)() != (*int_struct)){
    std::cout << "Fail 2" << std::endl;
  }

  return 0;
}
