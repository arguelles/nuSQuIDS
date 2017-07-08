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
#include "nuSQuIDS/nuSQuIDS.h"
#include "NSI.h"

/*
 * This example is exactly the same as the nuSQUIDSatm case but here 
 * we are using a different nuSQuIDS base class, in this concrete case 
 * we use the NSI nuSQuIDS subclass defined in the NSI example in the file
 * NSI.h
 */


using namespace nusquids;


//Function that gives the initial flux, for this example we set it to 1.0
double flux_function(double enu, double cz){
  return 1.;
}


int main()
{
  //Units and constants class
  squids::Const units;
  //Number of neutrinos
  unsigned int numneu = 3;


  //Minimum and maximum values for the energy and cosine zenith, notice that the energy needs to have the 
  //units, if they are omitted the input is in eV.
  double Emin=1.e0*units.GeV;
  double Emax=1.e3*units.GeV;
  double czmin=-1;
  double czmax=0;
  //Value of epsilon_mutau
  //double epsilon_mutau=1e-2;
  double epsilon_mutau=0.01;
  //Declaration of the atmospheric object now with the nuSQUIDSNSI instead of the basic nuSQuIDS
  std::cout << "Begin: constructing nuSQuIDS-Atm object" << std::endl;
  nuSQUIDSAtm<nuSQUIDSNSI> nus_atm(linspace(czmin,czmax,40),epsilon_mutau,logspace(Emin,Emax,600),numneu,both,false);
  std::cout << "End: constructing nuSQuIDS-Atm object" << std::endl;

  //Optionally you can access every one of the nuSQUIDSNSI objects in the array of zenith angles to 
  //for example set parameters, here is an example of a loop that sets the value of epsilon_mutau  
  for(nuSQUIDSNSI& nsq : nus_atm.GetnuSQuIDS()){
    nsq.Set_mutau(epsilon_mutau);
  }
  
  std::cout << "Begin: setting mixing angles." << std::endl;
  // set mixing angles, mass differences and cp phases
  nus_atm.Set_MixingAngle(0,1,0.59);
  nus_atm.Set_MixingAngle(0,2,0.154085);
  //nus_atm.Set_MixingAngle(1,2,0.6847);
  nus_atm.Set_MixingAngle(1,2,0.68479);
  
  nus_atm.Set_SquareMassDifference(1,7.65e-05);
  //nus_atm.Set_SquareMassDifference(1,7.54e-05);
  nus_atm.Set_SquareMassDifference(2,0.00243);
  
  nus_atm.Set_CPPhase(0,2,0);
  if(numneu > 3){
    nus_atm.Set_SquareMassDifference(3,1.);
    nus_atm.Set_MixingAngle(1,3,0.5);
  }
  std::cout << "End: setting mixing angles." << std::endl;
  
  // Setup integration precision
  nus_atm.Set_rel_error(1.0e-6);
  nus_atm.Set_abs_error(1.0e-6);
  nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);  

  //Array that contains the values of the energies and cos(zenith)
  auto e_range = nus_atm.GetERange();
  auto cz_range = nus_atm.GetCosthRange();
  
  std::cout << "Begin: setting initial state." << std::endl;
  //Construct the initial state, we set a flat spectra in zenith and log-energy
  marray<double,4> inistate{nus_atm.GetNumCos(),nus_atm.GetNumE(),2,numneu};
  std::fill(inistate.begin(),inistate.end(),0);
  for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
    for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
      for ( int rho = 0; rho < 2; rho ++ ){
        for (int flv = 0; flv < numneu; flv++){
          inistate[ci][ei][rho][flv] = (flv == 1) ? flux_function(e_range[ei], cz_range[ci]) : 0.0;//set 1 only to the muon flavor
        }
      }
    }
  }

  //Set the initial state in the atmSQuIDS object
  nus_atm.Set_initial_state(inistate,flavor);
  std::cout << "End: setting initial state." << std::endl;


  //Set to true the monitoring prgress bar and the vacuum oscillations
  nus_atm.Set_ProgressBar(true);
  nus_atm.Set_IncludeOscillations(true);


  //Here we evolve the system
  std::cout << "Begin: Evolution" << std::endl;
  nus_atm.EvolveState();
  std::cout << "End: Evolution" << std::endl;
  // We can save the current state in HDF5 format for future use.
  //nus_atm.WriteStateHDF5("state.hdf5");
  //For reading the hdf5 file that contains all the state, you can use 
  //nus_atm.ReadStateHDF5("state.hdf5")

  //This file will contine the final flux, since initially we set it to 1 for the muon, 
  //this can be readed as the muon ration fin/init in zenith and energy.
  std::ofstream file("fluxes_flavor.txt");


  //Set the resolution and the ranges for the output, remember that an interpolation in energy is used in 
  //in the interaction picture, the vacuum oscillations are solve analytically with arbitrary precision.
  int Nen =700;
  int Ncz=100;
  double lEmin=log10(Emin);
  double lEmax=log10(Emax);;

  //Writing to the file!  
  file << "# log10(E) cos(zenith) E flux_i . . . ." << std::endl;
  //for(double cz=czmin;cz<czmax;cz+=(czmax-czmin)/(double)Ncz){
  double cz = -1.;
    for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
      double E=pow(10.0,lE);
      file << lE << " " << cz << " " << E;
      for(int fl=0; fl<numneu; fl++){
	file << " " <<  nus_atm.EvalFlavor(fl,cz, E,0) << " " << nus_atm.EvalFlavor(fl,cz, E,1);
      }
      file << std::endl;
    }
    //file << std::endl;
  //}

  //Code to ask if you want to run the plotting script..
  std::string plt;
  std::cout << std::endl <<  "Done! " << std::endl <<
    "  Do you want to run the gnuplot script? yes/no" << std::endl;
  std::cin >> plt;
  if(plt=="yes" || plt=="y")
    return system("./plot.plt");

  return 0;
}
