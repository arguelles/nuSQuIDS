#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <chrono>  // for high_resolution_clock

using namespace nusquids;

int main(){

  squids::Const units;
//   unsigned int evalThreads = 8;
//   std::vector<nuSQUIDS> nusq_array;
  
  marray<double,2> lengths {8,10};
  marray<double,2> densities {8,10};
  marray<double,2> ye {8,10};
  marray<double,1> energies {8};
  // These lengths and densities are completely made up nonsense.
  // The intention is that we would have a pybinding in the real
  // world where we pass explicit distances and densities that were
  // calculated elsewhere.
  for (int i = 0; i < lengths.extent(0); i++){
    for (int j = 1; j < lengths.extent(1) + 1; j++){
      lengths[i][j-1] = 5000./j * units.km;
      densities[i][j-1] = (13. + i)/j;
      ye[i][j-1] = 0.5;
    }
    energies[i] = 10.*units.GeV;
  }
  
  std::cout << "Begin: constructing nuSQuIDS-Layers object" << std::endl;
  nuSQUIDSLayers<> nus_layer(lengths, densities, ye, energies, 3, neutrino);
  std::cout << "End: constructing nuSQuIDS-Layers object" << std::endl;
  
  nus_layer.Set_MixingParametersToDefault();
  nus_layer.Set_AllowConstantDensityOscillationOnlyEvolution(false);
  
  marray<double,1> ini_state({3},{0,1,0});
  nus_layer.Set_initial_state(ini_state, flavor);
  
  marray<double,2> states;
  
  std::cout << "Initial interaction picture states:" << std::endl;
  states = nus_layer.GetStatesArr();
  for (int i=0; i<states.extent(0); i++){
   for (int j=0; j<states.extent(1); j++)
     std::cout << states[i][j] << "   ";
   std::cout << std::endl;
  }
  
  std::cout << "Initial flavor state:" << std::endl;
  for (int n = 0; n < lengths.extent(0); n ++){
    std::cout << n << ": ";
    for(int i = 0; i < 3; i++){
      std::cout << nus_layer.EvalFlavorAtNode(i, n) << " ";
    }
    std::cout << std::endl;
  }
  
  nus_layer.Set_EvalThreads(1);
  
  auto start = std::chrono::high_resolution_clock::now();
  
  nus_layer.EvolveState();
  
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "evolution time: " << elapsed.count() << " s" << std::endl;
  
  std::cout << "Final interaction picture states:" << std::endl;
  states = nus_layer.GetStatesArr();
  for (int i=0; i<states.extent(0); i++){
   for (int j=0; j<states.extent(1); j++)
     std::cout << states[i][j] << "   ";
   std::cout << std::endl;
  }
  
  std::cout << "Final flavor state:" << std::endl;
  for (int n = 0; n < lengths.extent(0); n ++){
    std::cout << n << ": ";
    for(int i = 0; i < 3; i++){
      std::cout << nus_layer.EvalFlavorAtNode(i, n) << " ";
    }
    std::cout << std::endl;
  }
  
  // The following demonstrates evaluation with a given interaction picture state.
  // We just use the states as they are computed at the nodes, but in real life we
  // would do some sort of interpolation externally.
  std::cout << "Evaluating with states:" << std::endl;
  marray<double,1> summed_lengths = lengths.sum(1);
  
  for (int n = 0; n < lengths.extent(0); n ++){
    std::cout << n << ": ";
    for(int i = 0; i < 3; i++){
      // Is there no way to convert the slice more easily?
      marray<double,1> state {states.extent(1)};
      for (int j=0; j<state.extent(0); j++)
        state[j] = states[n][j];
      std::cout << nus_layer.EvalWithState(
        i, summed_lengths[n], 10.*units.GeV, state
      ) << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
