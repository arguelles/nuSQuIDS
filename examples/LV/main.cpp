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
#include "LV.h"


/*
 * This example it illustrates how to use the create a derived class for nuSQuIDS
 * in order to add some new physics, in this case Lorentz Violation (LV).
 * In the example we will create two objects  One with no LV and one with with non-zero value of c_mutau.
 * We calculate the muon neutrino disapperance for 10000 km and plot the results.
 */

using namespace nusquids;

int main()
{
  //class that contains the basic constant and units
  squids::Const units;
  //number of neutrinos
  unsigned int numneu=3;
  //Value for the epsilon mutau
  double re_c_mutau =1.0e-23*units.GeV;
  //minimum and maximum energy
  double Emin=1.e1*units.GeV;
  double Emax=1.e3*units.GeV;
  //Declaration of the nusquids LV object, for more details about how to construct an object
  //like this look at the file LV.h
  nuSQUIDSLV nus_lv(logspace(Emin,Emax,200),numneu,antineutrino,false);
  nuSQUIDSLV nus_std(logspace(Emin,Emax,200),numneu,antineutrino,false);

  gsl_complex zero {0.0*units.eV,0.0*units.eV};
  LVParameters null {zero, zero};
  nus_std.Set_LV_OpMatrix(null);
  nus_lv.Set_LV_EnergyPower(0.0);

  gsl_complex cmutau {re_c_mutau*units.eV,0.0*units.eV};
  LVParameters non_null_lv{ zero, cmutau };
  nus_lv.Set_LV_OpMatrix(non_null_lv);
  nus_lv.Set_LV_EnergyPower(0.0);

  //Setting up the propagation medium and track. For the track the number provided is the baseline.
  std::shared_ptr<Earth> earth = std::make_shared<Earth>();
  std::shared_ptr<Earth::Track> track = std::make_shared<Earth::Track>(10000.*units.km);

  nus_std.Set_Body(earth);
  nus_std.Set_Track(track);

  nus_lv.Set_Body(earth);
  nus_lv.Set_Track(track);

  //Set mixing angles and masses
  nus_std.Set_MixingAngle(0,1,0.563942);
  nus_std.Set_MixingAngle(0,2,0.154085);
  nus_std.Set_MixingAngle(1,2,0.785398);

  nus_std.Set_SquareMassDifference(1,7.65e-05);
  nus_std.Set_SquareMassDifference(2,0.00247);

  nus_lv.Set_MixingAngle(0,1,0.563942);
  nus_lv.Set_MixingAngle(0,2,0.154085);
  nus_lv.Set_MixingAngle(1,2,0.785398);

  nus_lv.Set_SquareMassDifference(1,7.65e-05);
  nus_lv.Set_SquareMassDifference(2,0.00247);

  //Setup integration settings
  nus_std.Set_rel_error(1.0e-15);
  nus_std.Set_abs_error(1.0e-15);

  nus_lv.Set_rel_error(1.0e-15);
  nus_lv.Set_abs_error(1.0e-15);

  marray<double,1> E_range = nus_std.GetERange();

  //Construct the initial state
  marray<double,2> inistate({200,numneu});
  double N0 = 1.0;
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        // initialize muon state
        inistate[i][k] = (k == 1) ? N0 : 0.0;
      }
  }

  // set the initial state
  nus_std.Set_initial_state(inistate,flavor);
  nus_lv.Set_initial_state(inistate,flavor);

  //Here we propagate the two neutrino fluxes, one with NSI and the other without NSI
  nus_lv.Set_ProgressBar(true);
  nus_std.Set_ProgressBar(true);
  std::cout <<"propagating the standard case, non-LV..." << std::endl;
  nus_std.EvolveState();
  std::cout << std::endl;

  std::cout <<"propagating the NSI case with c_mutau=" << re_c_mutau << " ..." << std::endl;
  nus_lv.EvolveState();
  std::cout << std::endl;

  int Nen =1000;
  double lEmin=log10(Emin);
  double lEmax=log10(Emax);

  std::ofstream file("fluxes_flavor.txt");

  file << "# log10(E) E flux_NSI_i flux_noNSI_i . . . ." << std::endl;
  for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
    double E=pow(10.0,lE);
    file << lE << " " << E << " ";
    for(int fl=0; fl<numneu; fl++){
      file << " " <<  nus_std.EvalFlavor(fl, E) << " " <<  nus_lv.EvalFlavor(fl, E);
    }
    file << std::endl;
  }
  std::string plt;
  std::cout << std::endl <<  "Done! " << std::endl <<
    "  Do you want to run the gnuplot script? yes/no" << std::endl;
  std::cin >> plt;
  if(plt=="yes" || plt=="y")
    return system("./plot.plt");

  return 0;
}
