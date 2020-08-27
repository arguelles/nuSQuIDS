#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <chrono>  // for high_resolution_clock

using namespace nusquids;

int main(){

  squids::Const units;
  unsigned int evalThreads = 1;
  std::vector<nuSQUIDS> nusq_array;
  
  marray<double,2> lengths {8,10};
  marray<double,2> densities {8,10};
  for (int i = 0; i < lengths.extent(0); i++){
    for (int j = 1; j < lengths.extent(1) + 1; j++){
      lengths[i][j-1] = 5000./j * units.km;
      densities[i][j-1] = (13. + i)/j;
    }
  }
  std::vector< std::vector< std::shared_ptr<ConstantDensity>>> const_dens_array;
  std::vector< std::vector< std::shared_ptr<ConstantDensity::Track>>> const_dens_track_array;
  
  for (int i = 0; i < lengths.extent(0); i++){
    std::vector< std::shared_ptr<ConstantDensity>> vdens;
    std::vector< std::shared_ptr<ConstantDensity::Track>> vtrk;
    
    for (int j = 0; j < lengths.extent(1); j++){
      vdens.push_back(std::make_shared<ConstantDensity>(densities[i][j], 0.5));
      vtrk.push_back(std::make_shared<ConstantDensity::Track>(lengths[i][j]));
    }
    const_dens_array.push_back(vdens);
    const_dens_track_array.push_back(vtrk);
  }
  
  marray<double,1> ini_state({3},{0,1,0});
  // construct nusquids objects
  for (int i = 0; i < lengths.extent(0); i++){
    nusq_array.emplace_back(3, neutrino);
    nusq_array.back().Set_E(10.*units.GeV);
    // set body and track to the first layer of the array of layers 
    // corresponding to this particular nusquids object.
    nusq_array.back().Set_Body(const_dens_array[i][0]);
    nusq_array.back().Set_Track(const_dens_track_array[i][0]);
    nusq_array.back().Set_MixingParametersToDefault();
    nusq_array.back().Set_AllowConstantDensityOscillationOnlyEvolution(true);
    nusq_array.back().Set_initial_state(ini_state, flavor);
  }

  std::cout << "Ini state in nus objects:" << std::endl;
  for (int n = 0; n < lengths.extent(0); n ++){
    std::cout << n << ": ";
    for(int i = 0; i < 3; i++){
      std::cout << nusq_array[n].EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  auto start = std::chrono::high_resolution_clock::now();
  
  if(evalThreads==1){
    std::cout << "using single task mode" << std::endl;
    for (int n = 0; n < lengths.extent(0); n++){
      std::cout << "evolving nusq object " << n << std::endl;
      for (int i = 0; i < lengths.extent(1); i++){
        nusq_array[n].Set_Body(const_dens_array[n][i]);
        nusq_array[n].Set_Track(const_dens_track_array[n][i]);
        nusq_array[n].EvolveState();
        //for(int j = 0; j < 3; j++){
        //  std::cout << nusq_array[n].EvalFlavor(j) << " ";
        //}
        //std::cout << std::endl;
      }
    }
  }
  else{
    std::cout << "using " << evalThreads << " threads" << std::endl;
    for (int i = 0; i < lengths.extent(1); i++){
      ThreadPool tpool(evalThreads);
      std::vector<std::future<void>> tasks;
      tasks.reserve(nusq_array.size());
      std::cout << "working on layer " << i << std::endl;
      for(int n = 0; n < nusq_array.size(); n++){
        nusq_array[n].Set_Body(const_dens_array[n][i]);
        nusq_array[n].Set_Track(const_dens_track_array[n][i]);
      }
      for(nuSQUIDS& nsq : nusq_array)
         tasks.emplace_back(tpool.enqueue([&](){ nsq.EvolveState(); }));
      for(const auto& task : tasks){
        task.wait();
      }
    }
    ////////////////////
    // It would probably be more efficient if every thread could do the entire 
    // evolution through all layers, rather than waiting for each layer to finish
    // in all nuSQuIDS calculators. I was trying to do it below, but it doesn't work...
    // If you have better C++ skills than me, you can probably make it work. 
    //
    // std::cout << "using " << evalThreads << " threads" << std::endl;
    // ThreadPool tpool(evalThreads);
    // std::vector<std::future<void>> tasks;
    // tasks.reserve(nusq_array.size());
    // 
    // for(int n = 0; n < nusq_array.size(); n++){
    //   nuSQUIDS& nus = nusq_array[n];
    //   std::vector< std::shared_ptr<ConstantDensity>> vdens = const_dens_array[n];
    //   std::vector< std::shared_ptr<ConstantDensity::Track>> vtrk = const_dens_track_array[n];
    //   
    //   tasks.emplace_back(tpool.enqueue([&](){
    //     for (int i = 0; i < vdens.size(); i++){
    //       nus.Set_Body(vdens[i]);
    //       nus.Set_Track(vtrk[i]);
    //       nus.EvolveState();
    //     }
    //   }));
    // }
    // for(const auto& task : tasks){
    //   task.wait();
    // }
    ///////////////////
  }

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "evolution time: " << elapsed.count() << " s" << std::endl;
  
  std::cout << "Final state in nus objects:" << std::endl;
  for (int n = 0; n < lengths.extent(0); n ++){
    std::cout << n << ": ";
    for(int i = 0; i < 3; i++){
      std::cout << nusq_array[n].EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
