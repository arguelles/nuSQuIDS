#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <time.h>

using namespace nusquids;

int main(){

  int numneu=3;
  squids::Const units;
  nusquids::marray<double,2> ND_flux_numode;
  nusquids::marray<double,2> ND_flux_anumode;
  nusquids::marray<double,2> FD_flux_numode;
  nusquids::marray<double,2> FD_flux_anumode;
  
  //reading the fluxes. Files are the reference beam from http://home.fnal.gov/~ljf26/DUNE2015CDRFluxes/ 
  ND_flux_numode=nusquids::quickread("./fluxes/g4lbne_v3r2p4b_FHC_ND_globes_flux.txt");
  ND_flux_anumode=nusquids::quickread("./fluxes/g4lbne_v3r2p4b_RHC_ND_globes_flux.txt");
  FD_flux_numode=nusquids::quickread("./fluxes/g4lbne_v3r2p4b_FHC_FD_globes_flux.txt");
  FD_flux_anumode=nusquids::quickread("./fluxes/g4lbne_v3r2p4b_RHC_FD_globes_flux.txt");
  
  //We set the nusquids object with the same grid in energy as the flux files, notice that this is not correct if the 
  //neutrino and antineutrino have different grid.
  int size=ND_flux_numode.size()/ND_flux_numode.extent(1);
  nuSQUIDS nus(linspace(ND_flux_numode[0][0]*units.GeV,ND_flux_numode[size-1][0]*units.GeV,size),numneu,both,false);
  nuSQUIDS anus(linspace(ND_flux_numode[0][0]*units.GeV,ND_flux_numode[size-1][0]*units.GeV,size),numneu,both,false);
  
  //Vector that contains all the energies.
  auto energy_range = nus.GetERange();

  //set mixing angles and masses to the default value.
  nus.Set_MixingParametersToDefault();

 
  // Set the density, electron fraction and the path where the neutirnos are going to travel
  const double lenght = 1297.*units.km;
  std::shared_ptr<ConstantDensity> matter = std::make_shared<ConstantDensity>(3.5,0.5); // density [gr/cm^3], ye [dimensionless]
  std::shared_ptr<ConstantDensity::Track> track = std::make_shared<ConstantDensity::Track>(lenght);

  // We set the configuration to nuSQuIDS
  nus.Set_Body(matter);
  nus.Set_Track(track);

  anus.Set_Body(matter);
  anus.Set_Track(track);
 
  
  // Construct the initial state
  marray<double,3> inistate_nus{nus.GetNumE(), 2, nus.GetNumNeu()};
  marray<double,3> inistate_anus{anus.GetNumE(), 2, anus.GetNumNeu()};

  for ( int ei = 0 ; ei < nus.GetNumE(); ei++){
    for ( int rho = 0; rho < 2; rho ++ ){
      for (int flv = 0; flv < numneu; flv++){
	inistate_nus[ei][rho][flv] = FD_flux_numode[ei][rho*3+flv+1];
	inistate_anus[ei][rho][flv] = FD_flux_anumode[ei][rho*3+flv+1];
      }
    }
  }

  // set initial state
  anus.Set_initial_state(inistate_anus,flavor);
  nus.Set_initial_state(inistate_nus,flavor);

  clock_t t;
  t=clock();
  // evolve the system through the first layer
  nus.EvolveState();
  anus.EvolveState();
  t=clock()-t;

  std::cout << "Propagation time: " << ((float)t)/CLOCKS_PER_SEC <<" seconds" << std::endl;


  //At this points both neutrino and antineutrino mode has the fluxes from the 
  //near detector propagated.

  //Print the fluxes, the columns are:
  //Energy FD_flux_numode_e FD_propflux_numode_e FD_flux_numode_mu FD_propflux_numode_mu
  //FD_flux_numode_tau FD_propflux_numode_tau FD_flux_numode_antie FD_propflux_numode_antie ....
  std::ofstream file("fluxes_flavor.txt");
  for(unsigned int i = 0 ; i < nus.GetNumE(); i++){
      double E = energy_range[i];
      file << E/units.GeV << " ";
      for(unsigned int k = 0; k < nus.GetNumNeu(); k ++){
	int rho=0;
        file << FD_flux_numode[i][rho*3+k+1] <<" " << nus.EvalFlavorAtNode(k,i,rho)
	     << " ";
      }
      for(unsigned int k = 0; k < nus.GetNumNeu(); k ++){
        int rho=1;
	file << FD_flux_numode[i][rho*3+k+1] <<" " << nus.EvalFlavorAtNode(k,i,rho)
	  << " ";
      }
      for(unsigned int k = 0; k < nus.GetNumNeu(); k ++){
	int rho=0;
        file << FD_flux_anumode[i][rho*3+k+1] <<" " << anus.EvalFlavorAtNode(k,i,rho)
	  << " ";
      }
      for(unsigned int k = 0; k < nus.GetNumNeu(); k ++){
	int rho=1;
        file << FD_flux_anumode[i][rho*3+k+1] <<" " << anus.EvalFlavorAtNode(k,i,rho)
	 << " ";
      }
      file << std::endl;
  }
  file.close();

  std::string plt;
  std::cout << std::endl <<  "Done! " << std::endl <<
    "  Do you want to run the gnuplot script? yes/no" << std::endl;
  std::cin >> plt;
  if(plt=="yes" || plt=="y"){
    return system("./plot.plt");
    std::cout << "The plots are in the ./plots/ folder" << std::endl;
  }
  return 0;
}
