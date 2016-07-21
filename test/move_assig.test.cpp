#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace nusquids;

double DC(double x){
  if (fabs(x) < 1.0e-15)
    return 0.;
  else
    return x;
}

int main(){
  squids::Const units;

  std::cout << std::setprecision(3);
  std::cout << std::fixed;

  nuSQUIDS nus(logspace(1.*units.GeV,1.e2*units.GeV,199),3,neutrino,false);

  std::shared_ptr<Earth> earth = std::make_shared<Earth>();
  std::shared_ptr<Earth::Track> track = std::make_shared<Earth::Track>(500.0*units.km);

  nus.Set_Body(earth);
  nus.Set_Track(track);

  // set mixing angles and masses
  nus.Set_MixingAngle(0,1,0.563942);
  nus.Set_MixingAngle(0,2,0.154085);
  nus.Set_MixingAngle(1,2,0.785398);
  nus.Set_SquareMassDifference(1,7.65e-05);
  nus.Set_SquareMassDifference(2,0.00247);

  // setup integration settings
  nus.Set_h_max( 500.0*units.km );
  nus.Set_rel_error(1.0e-9);
  nus.Set_abs_error(1.0e-9);

  marray<double,1> E_range = nus.GetERange();
  marray<double,2> inistate{E_range.size(),3};
  double N0 = 1.0e18;
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        inistate[i][k] = (k == 1) ? N0*pow(E_range[i],-2) : 0.0;
      }
  }
  // set the initial state
  nus.Set_initial_state(inistate,flavor);

  nuSQUIDS new_object(std::move(nus));

  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        std::cout << DC(fabs(inistate[i][k] - new_object.EvalFlavorAtNode(k,i))) << " ";
      }
      std::cout << std::endl;
  }

 new_object.EvolveState();
 marray<double,2> finalstate{E_range.size(),3};
 for ( int i = 0 ; i < finalstate.extent(0); i++){
   for ( int k = 0; k < finalstate.extent(1); k ++){
     finalstate[i][k] = new_object.EvalFlavorAtNode(k,i);
   }
 }
 nuSQUIDS other_object = std::move(new_object);

 for ( int i = 0 ; i < finalstate.extent(0); i++){
   for ( int k = 0; k < finalstate.extent(1); k ++){
     std::cout << fabs(finalstate[i][k] - other_object.EvalFlavorAtNode(k,i)) << " ";
   }
   std::cout << std::endl;
 }

 return 0;
}
