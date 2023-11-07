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

/*
 * The problem of solving the propagation of the atmospheric neutrinos is an
 * energy and zenith dependent problem, for this we include the class nuSQUIDSAtm
 * that allows to solve a set of nuSUIDS energy dependent objects to take in to account the 
 * zenith dependence. The problem of flux insertion into the atmosphere due to cosmic ray sources
 * during evolution requires a class derived from EarthAtm to be called into nuSQUIDSAtm
 */


using namespace nusquids;


//Function that sets the initial flux if you want to add flux in at the top of the atmosphere
//for this example we set it to 0
double flux_function(double enu, double cz){
  return 0;
}


int main()
{
  //Units and constants class
  squids::Const units;
  //Number of neutrinos (3) standard 4 1-sterile
  //Number of neutrinos, 4 is the case with one sterile neutrino
  std::cout << "Enter number of neutrino flavors to consider:\n";
  std::cout << "(3) Three Active Neutrinos, " << "(4) 3+1 Three Active and One Sterile Neutrino" << std::endl;
  unsigned int numneu;
  std::cin >>numneu;
  if( not(numneu==3 || numneu==4)){
    throw std::runtime_error("Only (3) or (4) are valid options");
  }
  std::cout << "Consider Earth absorption and interactions? " << std::endl;
  std::string use_int;
  std::cin >> use_int;
  bool interactions = false;
  if(use_int=="yes" || use_int=="y")
    interactions = true;

  //Minimum and maximum values for the energy and cosine zenith, notice that the energy needs to have the 
  //units, if they are omitted the input is in eV.
  double Emin=1e-2*units.GeV;
  double Emax=1.e6*units.GeV;
  double czmin=-1;
  double czmax=1;
  
  //Declaration of the atmospheric object
  //Must include the nuSQUIDS template or nuSQUIDSNSI user defined class derived from nuSQUIDS in order to 
  //use the EmittingEarthAtm class or another user defined class derived from EarthAtm
  std::cout << "Begin: constructing nuSQuIDS-Atm object" << std::endl;
  nuSQUIDSAtm<nuSQUIDS,EmittingEarthAtm> nus_atm(linspace(czmin,czmax,40),logspace(Emin,Emax,100),numneu,both,interactions);
  std::cout << "End: constructing nuSQuIDS-Atm object" << std::endl;
  
  //Sets height of atmosphere in km (defaults to 22 km when not set)
  nus_atm.Set_AtmHeight(60);

  //Functions Necessary for injection
  //Passes the energy range into the emitting bodies for interpolation
  nus_atm.Set_AtmEmissionEnergies();
  //Allows the nusquids objects to accept neutrino sources
  nus_atm.Set_NeutrinoSources(true);
  
  std::cout << "Begin: setting mixing angles." << std::endl;
  // set mixing angles, mass differences and cp phases
  nus_atm.Set_MixingAngle(0,1,0.563942);
  nus_atm.Set_MixingAngle(0,2,0.154085);
  nus_atm.Set_MixingAngle(1,2,0.785398);
  
  nus_atm.Set_SquareMassDifference(1,7.65e-05);
  nus_atm.Set_SquareMassDifference(2,0.00247);
  
  nus_atm.Set_CPPhase(0,2,0);
  if(numneu > 3){
    nus_atm.Set_SquareMassDifference(3,-1.);
    nus_atm.Set_MixingAngle(1,3,0.160875);
  }
  std::cout << "End: setting mixing angles." << std::endl;
  
  //Setup integration precision
  nus_atm.Set_rel_error(1.0e-6);
  nus_atm.Set_abs_error(1.0e-6);
  nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);  

  //Array that contains the values of the energies and cosine of the zenith, is the same length for every zenith
  auto e_range = nus_atm.GetERange();
  auto cz_range = nus_atm.GetCosthRange();
  
  std::cout << "Begin: setting initial state." << std::endl;

  //Construct the initial state at the top of the atmosphere, we set it to 0 here
  marray<double,4> inistate{nus_atm.GetNumCos(),nus_atm.GetNumE(),2,numneu};
  std::fill(inistate.begin(),inistate.end(),0);
  for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
    for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
      for ( int rho = 0; rho < 2; rho ++ ){
        for (int flv = 0; flv < numneu; flv++){
          inistate[ci][ei][rho][flv] = (flv == 0) ? flux_function(e_range[ei], cz_range[ci]) : 0;//set the flux for the muon flavor: 1
          inistate[ci][ei][rho][flv] = (flv == 1) ? flux_function(e_range[ei], cz_range[ci]) : 0;//set the flux for the muon flavor: 1
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

  //Here we do the evolution of all the states
  std::cout << "Begin: Evolution" << std::endl;
  nus_atm.EvolveState();
  std::cout << "End: Evolution" << std::endl;

  //We can save the current state in HDF5 format for future use.
  //nus_atm.WriteStateHDF5("./atmospheric_emission_example_numneu_"+std::to_string(numneu)+".hdf5");

  //This file will contain the final flux, since initially we set it to 1 for the muon, 
  //this can be read as the muon ration F_final/F_initial in cos(zenith) and energy.
  std::ofstream file("fluxes_flavor.txt");


  //Set the resolution and the ranges for the ouput, remember that an interpolation in energy is used in 
  //in the interaction picture, the vacuum oscillations are solve analytically with arbitrary Energy precision.
  //For the zenith a linear interpolation is used.
  int Nen=700;
  int Ncz=100;
  double lEmin=log10(Emin);
  double lEmax=log10(Emax);;

  //Writing to the file!  
  file << "# log10(E) cos(zenith) E flux_i . . . ." << std::endl;
  for(double cz=czmin;cz<czmax;cz+=(czmax-czmin)/(double)Ncz){
    for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
      double E=pow(10.0,lE);
      file << lE - log10(units.GeV) << " " << cz << " " << E/units.GeV;
      for(int fl=0; fl<numneu; fl++){
        file << " " <<  nus_atm.EvalFlavor(fl,cz, E, 0);
      }
      for(int fl=0; fl<numneu; fl++){
        file << " " <<  nus_atm.EvalFlavor(fl,cz, E, 1);
      }
      file << std::endl;
    }
    file << std::endl;
  }

  //This asks if you want to run the gnuplot plotting script.
  std::string plt;
  std::cout << std::endl <<  "Done! " << std::endl <<
    "  Do you want to run the gnuplot script? yes/no" << std::endl;
  std::cin >> plt;
  if(plt=="yes" || plt=="y")
    return system("./plot.plt");


  return 0;
}
