#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace nusquids;

int main(){

  squids::Const units;

  nuSQUIDS nus(1.e0*units.GeV,1.e1*units.GeV,60,3,neutrino,false,false);
  auto energy_range = nus.GetERange();

  // set mixing angles and masses
  nus.Set_MixingParametersToDefault();

  // first layer
  const double layer_1 = 100.*units.km;
  std::shared_ptr<Vacuum> vacuum = std::make_shared<Vacuum>();
  std::shared_ptr<Vacuum::Track> track_env0 = std::make_shared<Vacuum::Track>(layer_1);
  // second layer
  const double layer_2 = 50.*units.km;
  std::shared_ptr<ConstantDensity> constdens_env1 = std::make_shared<ConstantDensity>(3.5,0.5); // density [gr/cm^3[, ye [dimensionless]
  std::shared_ptr<ConstantDensity::Track> track_env1 = std::make_shared<ConstantDensity::Track>(layer_2);
  // three layer
  const double layer_3 = 200.*units.km;
  std::shared_ptr<ConstantDensity> constdens_env2 = std::make_shared<ConstantDensity>(10.,0.1); // density [gr/cm^3[, ye [dimensionless]
  std::shared_ptr<ConstantDensity::Track> track_env2 = std::make_shared<ConstantDensity::Track>(layer_3);

  // set the first layer
  nus.Set_Body(vacuum);
  nus.Set_Track(track_env0);

  // construct the initial state
  marray<double,2> inistate{nus.GetNumE(),nus.GetNumNeu()};
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        inistate[i][k] = (k == 1);
      }
  }

  // set initial state
  nus.Set_initial_state(inistate,flavor);

  // evolve the system through the first layer
  nus.EvolveState();

  // go through next regions
  nus.Set_Body(constdens_env1);
  nus.Set_Track(track_env1);
  nus.EvolveState();

  nus.Set_Body(constdens_env2);
  nus.Set_Track(track_env2);
  nus.EvolveState();

  // print the oscillation probabilities
  std::ofstream file("fluxes_flavor.txt");
  for(unsigned int i = 0 ; i < nus.GetNumE(); i++){
      double E = energy_range[i];
      file << E/units.GeV << " ";
      for(unsigned int k = 0; k < nus.GetNumNeu(); k ++){
        double p = nus.EvalFlavorAtNode(k,i);
        file << p << " ";
      }
      file << std::endl;
  }
  file.close();

  std::string plt;
  std::cout << std::endl <<  "Done! " << std::endl <<
    "  Do you want to run the gnuplot script? yes/no" << std::endl;
  std::cin >> plt;
  if(plt=="yes" || plt=="y")
    return system("./plot.plt");

  return 0;
}
