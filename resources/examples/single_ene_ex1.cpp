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
 ******************************************************************************/

#include <vector>
#include <iostream>
#include "nuSQUIDS.h"

/*
 * This file demonstrates how to calculate neutrino oscillation
 * probabilities for different initial flavor states, mixing angles
 * and bodies.
 */

using namespace nusquids;

int main()
{
  // We must first create a nuSQUIDS object. In order to do this
  // we must specify the number of neutrino flavors and if we are
  // going to consider neutrino or antineutrino oscillations.

  // In this example we set N_neutrino = 3 and Type = "neutrino".
  nuSQUIDS nus(3,"neutrino");

  // We will now set the mixing paremeters and square mass
  // differences. The angles will be given in radians and the
  // differences in eV^2.

  // We use the standard parametrization as described in the
  // documentation.

  // mixing angles
  nus.Set(TH12,0.563942);
  //nus.Set("th12",0.);
  nus.Set(TH13,0.154085);
  //nus.Set(TH13,0.);
  nus.Set(TH23,0.785398);
  //nus.Set(TH23,0.);
  // square mass differences
  nus.Set(DM21SQ,7.65e-05);
  nus.Set(DM31SQ,0.00247);
  // CP phase
  nus.Set(DELTA1,0.0);

  // Now we set the neutrino energy which we are interested on.
  // Energies are always given in natural units. To handle
  // the units the nuSQUIDS object has a unit subclass
  // which contains the most used units.

  nus.Set_E(10.0*nus.units.GeV);

  // Next we have to specify where the neutrino oscillation
  // will take place and its trayectory. We call this, respectively,
  // body and track. nuSQUIDS comes with a set of predefined objects
  // which its specified trayectories. We will demonstrate its use
  // in this program.

  /*
   * Example 1
   * =========
   * A long baseline oscillation probability
   */
  std::cout << "*** Earth LongBaseline Neutrino Osc ***" << std::endl;

  double baseline = 500.0*nus.units.km;
  // create a body object, in this case the Earth
  std::shared_ptr<Earth> earth = std::make_shared<Earth>();
  // create a trajectory on the body. In this case the Earth
  // trajectory is given by three quantities: Initial position,
  // Final position, and Baseline.
  std::shared_ptr<Earth::Track> earth_track = std::make_shared<Earth::Track>(0.0,baseline,baseline);

  // Then we set the body and trajectory in the nuSQUIDS object.
  nus.Set_Body(earth);
  nus.Set_Track(earth_track);

  // Now we have to set the initial neutrino state. The state
  // is specified by a vector of doubles the size of this vector
  // must be the number of neutrino states. nuSQUIDS can initialize
  // the problem in either the mass or flavor basis. Thus, e.g.
  //
  // (1,0,0)
  //
  // corresponds to the nu_e flavor or nu_1 depending if the flavor
  // or mass basis is choosen. Similarly,
  //
  // (0,1,0) -> nu_mu (nu_2)
  // (0,0,1) -> nu_tau (nu_3)
  //
  // if more flavors are added then those are considered nu_sterile
  // and by default have no coherent matter interactions. Clearly,
  // one can also start with superposition of flavor or mass states
  // e.g. (1,1,1), (1,1,0), etc.

  // In this case we will start with a pure nu_mu state.
  std::vector<double> ini_state{0,1,0};
  nus.Set_initial_state(ini_state,"flavor");

  // Lets print out the initial state
  std::cout << "In state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  // Next we can set the numerical accuracy of our result. The values
  // of this parameters depend on the problem and hand, we encourage
  // the user to try different values.
  nus.Set_h_max( 200.0*nus.units.km );
  nus.Set_rel_error(1.0e-12);
  nus.Set_abs_error(1.0e-12);

  // Now that we have all the pieces in place we can tell the
  // nuSQUIDS object to evolve the given state.
  nus.EvolveState();

  // Finally, lets output the state result and save it.
  std::cout << "Out state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  /*
   * Example 2
   * =========
   * Atmospheric neutrino oscillation probability
   */
  double phi = acos(-1.0);

  // To calculate atmospheric neutrino oscillation probabilities
  // we need to speciefy a different body.

  std::cout << "*** Earth Atmospheric Neutrino Osc ***" << std::endl;

  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> earth_atm_track = std::make_shared<EarthAtm::Track>(phi);

  nus.Set_Body(earth_atm);
  nus.Set_Track(earth_atm_track);

  // We reset the initial condition
  nus.Set_initial_state(ini_state,"flavor");
  // We can change the energy
  nus.Set_E(5.0*nus.units.GeV);

  std::cout << "In state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  // We do the calculation
  nus.EvolveState();

  // Output the result
  std::cout << "Out state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  /*
   * Example 3
   * =========
   * Vacuum neutrino oscillation probability
   */

  std::cout << "*** Vacuum Neutrino Osc ***" << std::endl;
  std::shared_ptr<Vacuum> vacuum = std::make_shared<Vacuum>();
  double baseline_2 = 500.0*nus.units.km;
  std::shared_ptr<Vacuum::Track> track_vac = std::make_shared<Vacuum::Track>(0.0,baseline_2);

  nus.Set_Body(vacuum);
  nus.Set_Track(track_vac);

  // Lets set the initial state to electron
  ini_state = {1,0,0};
  nus.Set_initial_state(ini_state,"flavor");
  // We can change the energy some MeV
  nus.Set_E(15.0*nus.units.MeV);

  std::cout << "In state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  nus.EvolveState();

  // Output the result
  std::cout << "Out state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  /*
   * Example 4
   * =========
   * Constant density neutrino oscillation probability
   */

  std::cout << "*** Constant Density Neutrino Osc ***" << std::endl;
  // We will define a constant density environment
  // which depends on the density (in gr/cm^3) and electron fraction
  double density = 100.0; // gr/cm^3
  double ye = 0.3;
  std::shared_ptr<ConstantDensity> constdens = std::make_shared<ConstantDensity>(density,ye);
  double baseline_3 = 500.0*nus.units.km;
  std::shared_ptr<ConstantDensity::Track> track_constdens = std::make_shared<ConstantDensity::Track>(0.0,baseline_3);

  nus.Set_Body(constdens);
  nus.Set_Track(track_constdens);

  // Another thing we can do is change the oscillation parameters
  // lets try the following

  nus.Set(TH13,0.5);
  nus.Set(TH23,0.5);
  nus.Set(DM31SQ,0.0247);

  // Lets set the initial state to electron. Note that we have
  // to set the mixing parameters *before* we set the initial state
  // in the flavor basis.
  ini_state = {1,0,0};
  nus.Set_initial_state(ini_state,"flavor");
  // We can change the energy some MeV
  nus.Set_E(10.0*nus.units.MeV);

  std::cout << "In state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  nus.EvolveState();

  // Output the result
  std::cout << "Out state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  /*
   * Example 5
   * =========
   * Variable density neutrino oscillation probability
   */

  std::cout << "*** Variable Density Neutrino Osc ***" << std::endl;
  // We will define a variable density environment
  // which depends on the density (in gr/cm^3) and electron fraction
  // at each point of the neutrino trayectory. The given lists
  // will be interpolated using a spline.

  std::vector<double> x_arr(40);
  std::vector<double> density_arr(40);
  std::vector<double> ye_arr(40);

  double size = 1000.0*nus.units.km;
  for(int i = 0; i < 40; i++){
    x_arr[i] = size*(i/40.);
    density_arr[i] = fabs(cos((double)i));
    ye_arr[i] = fabs(sin((double)i));
  }

  std::shared_ptr<VariableDensity> vardens = std::make_shared<VariableDensity>(x_arr,density_arr,ye_arr);
  std::shared_ptr<VariableDensity::Track> track_vardens = std::make_shared<VariableDensity::Track>(0.0,200.0*nus.units.km);

  nus.Set_Body(vardens);
  nus.Set_Track(track_vardens);

  // Lets set the initial state to electron
  ini_state = {0,0,1};
  nus.Set_initial_state(ini_state,"flavor");
  // We can change the energy some MeV
  nus.Set_E(100.0*nus.units.GeV);

  std::cout << "In state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  nus.EvolveState();

  // Output the result
  std::cout << "Out state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  return 0;

  /*
   * Example 6
   * =========
   * Solar neutrino oscillation probability
   */

  std::cout << "*** Solar Neutrino Osc ***" << std::endl;
  // We will calculate the solar neutrino oscillation probability
  // using the standard solar model

  std::shared_ptr<Sun> sun = std::make_shared<Sun>();
  std::shared_ptr<Sun::Track> track_sun = std::make_shared<Sun::Track>(0.0,sun->radius*nus.units.km);

  nus.Set_Body(sun);
  nus.Set_Track(track_sun);

  // Lets set the initial state to electron
  ini_state = {1,0,0};
  nus.Set_initial_state(ini_state,"flavor");
  // We can change the energy some MeV
  nus.Set_E(10.0*nus.units.MeV);

  std::cout << "In state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }

  nus.EvolveState();

  // Output the result
  std::cout << "Out state" << std::endl;
  for (double EE : nus.GetERange()){
    std::cout << EE/nus.units.GeV << " ";
    for(int i = 0; i < 3; i++){
      std::cout << nus.EvalFlavor(i) << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
