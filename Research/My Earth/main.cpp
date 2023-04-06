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
  return 1.0;
}



int main()
{
  //Units and constants class
  squids::Const units;
  //Number of neutrinos
  unsigned int numneu = 7;


//Minimum and maximum values
  double Emin=1.e0*units.GeV;
  double Emax=1.e6*units.GeV;
  double czmin=-1;
  double czmax=0;

  //Declaration of the atmospheric object
  std::cout << "Begin: constructing nuSQuIDS-Atm object" << std::endl;
  nuSQUIDSAtm<nuSQUIDSNSI> nus_atm(linspace(czmin,czmax,40),logspace(Emin,Emax,100),numneu,both,false);
  std::cout << "End: constructing nuSQuIDS-Atm object" << std::endl;


  //1 or 0 to turn off various parts of model
  double sterile = 1;
  double quasi = 0;
  double standard = 1;
  
  double DMDif[numneu] = {quasi*QM0, quasi*QM1, quasi*QM2, sterile*SM3, standard*EM4, standard*MM5, standard*TM6};

  //mixing angles are written as QS1-E, QS2-M, QS3-T, VS-E, E-M, E-T, M-T
  double DMAng[numneu] = {quasi*Th04, quasi*Th15, quasi*Th26, sterile*Th34, standard*Th45, standard*Th46, standard*Th56};

  for(int i=1; i<numneu; i++){
  nus_atm.Set_SquareMassDifference(i,DMDif[i]-DMDif[0]);
  }
  //set mixing angles, mass differences and cp phases
  std::cout << "Begin: setting mixing angles." << std::endl;

  nus_atm.Set_MixingAngle(0,4,DMAng[0]);
  nus_atm.Set_MixingAngle(1,5,DMAng[1]);
  nus_atm.Set_MixingAngle(2,6,DMAng[2]);
  nus_atm.Set_MixingAngle(3,4,DMAng[3]);
  nus_atm.Set_MixingAngle(4,5,DMAng[4]);
  nus_atm.Set_MixingAngle(4,6,DMAng[5]);
  nus_atm.Set_MixingAngle(5,6,DMAng[6]);
  std::cout << "End: setting mixing angles." << std::endl;
  

 
  // Setup integration precision
  nus_atm.Set_rel_error(1.0e-6*100);
  nus_atm.Set_abs_error(1.0e-6*100);
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
          inistate[ci][ei][rho][flv] = (flv == 5) ? flux_function(e_range[ei], cz_range[ci]) : 0.0;//set 1 only to the muon flavor
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

  //This file will contine the final flux, since initially we set it to 1 for the muon, 
  //this can be readed as the muon ration fin/init in zenith and energy.
  std::ofstream file("fluxes_flavor.txt");


  //Set the resolution and the ranges for the ouput, remember that an interpolation in energy is used in 
  //in the interaction picture, the vacuum oscillations are solve analytically with arbitrary Energy precision.
  //For the zenith a linear interpolation is used.
  int Nen=700;
  int Ncz=100;
  double lEmin=log10(Emin);
  double lEmax=log10(Emax);

  //Writing to the file!  
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


  //This ask if you want to run the gnuplot plotting script.
  std::string plt;
  std::cout << std::endl <<  "Done! " << std::endl <<
    "  Do you want to run the gnuplot script? yes/no" << std::endl;
  std::cin >> plt;
  if(plt=="yes" || plt=="y")
    return system("./plot.plt");

  return 0;
}
