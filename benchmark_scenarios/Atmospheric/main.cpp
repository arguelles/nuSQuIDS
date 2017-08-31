 /******************************************************************************
 *    This program is free software: you can redistribute it and/or modify     *
 *   it under the terms of the GNU General Public License as published by      *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   This program is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                             *
 *   Authors:                                                                  *
 *      Carlos Arguelles (University of Wisconsin Madison)                     *
 *         carguelles@icecube.wisc.edu                                         *
 *      Jordi Salvado (University of Wisconsin Madison)                        *
 *         jsalvado@icecube.wisc.edu                                           *
 *      Christopher Weaver (University of Wisconsin Madison)                   *
 *         chris.weaver@icecube.wisc.edu                                       *
 ******************************************************************************/

#include <vector>
#include <iostream>
#include <fstream>
#include <nuSQuIDS/nuSQuIDS.h>
#include <time.h>

/*
 * This example illustrates how to use the atmospheric example in a more practical case.
 * We implement the case of low atmospheric neutrinos where the scattering interactions 
 * can be neglected.
 */


using namespace nusquids;

//Class that reads the initial fluxes from a file in the folder ./fluxes
//The constructor needs to be called with the name of the model
class fluxes{
  nusquids::marray<double,2> flux;
  int size;
  int sizeZ,sizeE;
  std::string fluxpath=std::string("./fluxes/");
public:
  double minzenith,maxzenith,minE,maxE;
  fluxes(std::string name){
    SetFlux(name);
  }
  void SetFlux(std::string name){
    flux=nusquids::quickread(fluxpath+name);
    size=flux.size()/flux.extent(1);
    minzenith=flux[0][0];
    maxzenith=flux[size-1][0];
    minE=flux[0][1];
    maxE=flux[size-1][1];
    for(int i=0;i<size;i++){
      if(flux[0][0]!=flux[i][0]){
	sizeE=i;
	sizeZ=size/sizeE;
	break;
      }
    }
  }

  //returns the value of the flux, a bit too simple but fast.
  double get_flux(double z,double E,unsigned int flv, unsigned int type){
    int iz=(int)((sizeZ-1)*((z-minzenith)/(maxzenith-minzenith)));
    int iE=(int)((sizeE-1)*((log10(E)-log10(minE))/(log10(maxE)-log10(minE))));
    if(iz<0 or iz>=sizeZ)
      std::cerr << "Zenith out of range: " << iz << "  size:" << sizeZ << "  val: " << z <<   std::endl;
    if(iE<0 or iE>=sizeE)
      std::cerr << "Energy out of range: " << iE << "  size:" << sizeE << "  val: " << E <<std::endl;
    int idx=100;
    if(flv==0)
      idx=5;
    else if(flv==1)
      idx=3;
    else
      return 0.0;
	
    return flux[iz*sizeE+iE][idx + type];
  }

};


struct nuisance{
  double delta_gamma;
  double norm;
  double nuanu_ratio;
  double gamma_pivotE=1e12;
  
};

void propagate_fluxes(std::string initial_flux_name, double dm2, double sqth,std::string& name)
{
  //Units and constants class
  squids::Const units;
  //We create a flux instance with for the flux name as an argument
  fluxes fl(initial_flux_name);

  //We keep the number of neutrinos to 3 in this case, but is easy to generalize.
  unsigned int numneu=3;
  //Absortion (NC and CC scattering interactions can be turn on and off here 
  bool interactions = false;

  //Minimum and maximum values for the energy and cosine zenith, notice that the energy needs to have the 
  //units, if they are omitted the input is in eV.
  double Emin=fl.minE*units.GeV;
  double Emax=fl.maxE*units.GeV;
  double czmin=-1;//fl.minzenith;
  double czmax=1;//fl.maxzenith;

  std::cout << czmin << " " << czmax << std::endl;

  //Declaration of the atmospheric object, since there is an uncertatny in the pion/kaon 
  //ration we need two fluxes one for kaon and pion fluxes
  std::cout << "Begin: constructing nuSQuIDS-Atm object" << std::endl;
  nuSQUIDSAtm<> nus_atm(linspace(czmin,czmax,20),logspace(Emin,Emax,100),numneu,both,interactions);
  std::cout << "End: constructing nuSQuIDS-Atm object" << std::endl;
  
  
  std::cout << "Begin: setting mixing angles." << std::endl;
  //set mixing angles, mass differences and cp phases
  nus_atm.Set_MixingAngle(0,1,0.563942);
  nus_atm.Set_MixingAngle(0,2,0.154085);
  //  nus_atm_pion.Set_MixingAngle(1,2,0.785398);
  nus_atm.Set_MixingAngle(1,2,asin(sqrt(sqth)));
  nus_atm.Set_SquareMassDifference(1,7.65e-05);
  //nus_atm_pion.Set_SquareMassDifference(2,0.00247);
  nus_atm.Set_SquareMassDifference(2,sqrt(dm2));
  nus_atm.Set_CPPhase(0,2,0);

  nus_atm.Set_rel_error(1.0e-4);
  nus_atm.Set_abs_error(1.0e-4);
  nus_atm.Set_GSL_step(gsl_odeiv2_step_rkf45);  

  //Array that contains the values of the energies and cosine of the zenith, is the same length for every zenith
  auto e_range = nus_atm.GetERange();
  auto cz_range = nus_atm.GetCosthRange();
  
  //Construct the initial state, we set a flat spectra in zenith and log-energy
  marray<double,4> inistate{nus_atm.GetNumCos(),nus_atm.GetNumE(),2,numneu};

  std::fill(inistate.begin(),inistate.end(),0);
  for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
    for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
      for ( int rho = 0; rho < 2; rho ++ ){
        for (int flv = 0; flv < numneu; flv++){
          inistate[ci][ei][rho][flv] = fl.get_flux(cz_range[ci], e_range[ei]/units.GeV,flv, rho);
	  //set the flux for muon and electron neutrinos and antineutrinos
        }
      }
    }
  }
  //Set the initial state in the atmSQuIDS object
  nus_atm.Set_initial_state(inistate,flavor);


  //Set to true the monitoring prgress bar and the vacuum oscillations
  nus_atm.Set_ProgressBar(false);
  nus_atm.Set_IncludeOscillations(true);
  std::cout << "End: setting initial state." << std::endl;

  //Here we do the evolution of all the states
  std::cout << "Begin: Evolution" << std::endl;
  nus_atm.EvolveState();
  std::cout << "End: Evolution" << std::endl;

  //We save the propagated states in the folder prop_fluxes, notcie that we are saving the full 
  //density matrix therefore any expectation value can be computed loading this file.

  name="./prop_fluxes/atmospheric_dmq"+std::to_string(dm2)+"_sqth"+std::to_string(sqth)+".hdf5";

  nus_atm.WriteStateHDF5(name);
}

//Function that returns the flux at the detector for a given cosin(zenith), energy and flavor.
//The extra argument are the pion and kaon fluxes, precomputed for example with the function propagate_fluxes
//and a structure with the nuisance parameter.
double flux_at_detector(double cz, double E, int fl,  nuSQUIDSAtm<>& nus_atm, nuisance &sys){
  double nu=nus_atm.EvalFlavor(fl,cz, E, 0);
  double anu=nus_atm.EvalFlavor(fl,cz, E, 1);
  
  return (nu+anu*sys.nuanu_ratio)*sys.norm*pow(E/sys.gamma_pivotE,sys.delta_gamma);

}

int main(){
  clock_t t;
  squids::Const units;
  //value of the physical parameters physical paraemter.
  //For a paramter scan this can be in a loop to compute all the cases that we need
  double dm2=2.2e-3;
  double sq2th=0.8;
  std::string name;
  
  //Propagation of the flux for the given paramters, the fluxes are going to be saved
  //in the folder prop_fluxes to be used later.
  //The first argument is the name of the flux model, the second and third the value of the 
  //physical parameters.

  t=clock();
  propagate_fluxes("tflux.txt", dm2, sq2th, name);
  t=clock()-t;

  std::cout <<std::endl << "************ times ************" << std::endl;
  std::cout << "Propagation time: " << ((float)t)/CLOCKS_PER_SEC <<" seconds" << std::endl;

  //This part usually it would be in a different main file, for a better performance is better to propagate 
  //the fluxes and then load them and recover the physical flux at the detector that will depend on 
  //a set of systematic nuisance parameters.
  //first we recover the fluxes. 
  nuSQUIDSAtm<> nus_atm(name); 

  //We create a nuisance parmeter data structure.
  nuisance sys;
  //and we set the values.
  sys.delta_gamma=0.0;
  sys.norm=1.;
  sys.nuanu_ratio=1.;


  //finally we can call the flux_at_detector function and get the value of the flux for the paramters.
  //this funtion is fast to evaluate and it should be used as the imput for a given experiment, 
  //weighting a MC set or convoluted with a given detector responce function. Here we just 
  //evaluate a single energy and zenith.
  double Energy=1.e1*units.GeV;
  double CosZenith=-0.7;
  int fl=1;
  
  t=clock();
  double flux=flux_at_detector(CosZenith, Energy, flavor, nus_atm,sys); 
  t=clock()-t;
  std::cout << "Evaluation time: " << ((float)t)/CLOCKS_PER_SEC << " seconds" << std::endl;
  std::cout << "Message to take home:"<< std::endl << " Evaluate more and propagate less!" <<std::endl;
  std::cout << "*******************************" << std::endl<<std::endl;
  //Is important to notice that as it is done the systematic parameters enter at evaluation 
  //any systematic error study can be done for a set of already propagated fluxes. 

  std::cout << "The flux for Cos(Zenith)="<< CosZenith <<  "   E=" << Energy/units.GeV 
	    << "GeV   flavor=" <<fl << " is: " 
	    << flux
	    << " carlitos units" << std::endl; 

}
