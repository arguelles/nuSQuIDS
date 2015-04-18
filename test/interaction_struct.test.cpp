#include <nuSQuIDS/nuSQUIDS.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace nusquids;

int main(){

  nuSQUIDS nus1(1.e2,1.e6,60,3,neutrino,true,true);
  auto e_range = nus1.GetERange();
  auto int_struct = nus1.GetInteractionStructure();

  marray<double,1> costh_range {{10},{-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,}};
  nuSQUIDSAtm<> nusatm(costh_range,e_range,3,neutrino,true,nullptr);

  marray<double,3> sigmacc_1 = nusatm.GetInteractionStructure()->sigma_CC;
  marray<double,3> sigmacc_2 = int_struct->sigma_CC;

  if(nusatm.GetInteractionStructure() != int_struct){
    std::cout << "Fail" << std::endl;
  }

  return 0;
}
