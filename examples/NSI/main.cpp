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
#include <nuSQuIDS/nuSQuIDS.h>
#include "NSI.h"


/*
 * This example it illustrates how to use the create a derived class for nuSQuIDS
 * in order to add some new physics, in this case Non Standard Interactions(NSI).
 * In the example we will create 2 object on with non cerlo value of epsilon_mutau 
 * and the other with zero, finally we plot the results.
 */

using namespace nusquids;

int main()
{
  //class that contains the basic constant and units
  squids::Const units;
  //number of neutrinos
  int numneu=3;
  //Value for the epsilon mutau
  double eps_mutau=1.0e-2;
  //minimum and maximum energy
  double Emin=1.e1*units.GeV;
  double Emax=1.e3*units.GeV;
  //Declaration of the nusquids NSI object, for more details about how to construct an object
  //like this look at the file NSI.h
  //We declare two objects "nus" with NSI value of eps_mutau and nus_zero with 0.0
  nuSQUIDSNSI nus(eps_mutau,logspace(Emin,Emax,200),numneu,antineutrino,false);
  nuSQUIDSNSI nus_zero(0.0,logspace(Emin,Emax,200),numneu,antineutrino,false);

  //zenith angle for which we propagate the neutrino flux, they go though the earth.
  double phi = acos(-1.);
  //Setting up the object, and the track, the second depend on the zenith angle value
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track_atm = std::make_shared<EarthAtm::Track>(phi);
  nus.Set_Body(earth_atm);
  nus.Set_Track(track_atm);
  
  nus_zero.Set_Body(earth_atm);
  nus_zero.Set_Track(track_atm);


  //Set mixing angles and masses
  nus.Set_MixingAngle(0,1,0.563942);
  nus.Set_MixingAngle(0,2,0.154085);
  nus.Set_MixingAngle(1,2,0.785398);

  nus.Set_SquareMassDifference(1,7.65e-05);
  nus.Set_SquareMassDifference(2,0.00247);

  nus_zero.Set_MixingAngle(0,1,0.563942);
  nus_zero.Set_MixingAngle(0,2,0.154085);
  nus_zero.Set_MixingAngle(1,2,0.785398);

  nus_zero.Set_SquareMassDifference(1,7.65e-05);
  nus_zero.Set_SquareMassDifference(2,0.00247);


  //Setup integration settings
  nus.Set_h_max( 200.0*units.km );
  nus.Set_rel_error(1.0e-15);
  nus.Set_abs_error(1.0e-15);

  nus_zero.Set_h_max( 200.0*units.km );
  nus_zero.Set_rel_error(1.0e-15);
  nus_zero.Set_abs_error(1.0e-15);


  marray<double,1> E_range = nus.GetERange();

  //Construct the initial state
  marray<double,2> inistate({200,3});
  double N0 = 1.0;
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        // initialize muon state
        inistate[i][k] = (k == 1) ? N0 : 0.0;
      }
  }

  // set the initial state
  nus.Set_initial_state(inistate,flavor);
  nus_zero.Set_initial_state(inistate,flavor);


  //Here we propagate the two neutrino fluxes, one with NSI and the other without NSI
  nus.Set_ProgressBar(true);
  nus_zero.Set_ProgressBar(true);
  std::cout <<"propagating the standard case, non-NSI ..." << std::endl;
  nus.EvolveState();
  std::cout << std::endl;

  std::cout <<"propagating the NSI case with epsilon_mutau=" << eps_mutau << " ..." << std::endl;
  nus_zero.EvolveState();
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
      file << " " <<  nus.EvalFlavor(fl, E) << " " <<  nus_zero.EvalFlavor(fl, E);
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
