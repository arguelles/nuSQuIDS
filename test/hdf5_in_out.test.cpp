#include <nuSQuIDS/nuSQUIDS.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace nusquids;

int main(){

  nuSQUIDS nus(1.e2,1.e6,60,3,neutrino,true,true);

  double phi = acos(-0.5);
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track_atm = std::make_shared<EarthAtm::Track>(phi);

  nus.Set_Body(earth_atm);
  nus.Set_Track(track_atm);

  // set mixing angles and masses
  nus.Set_MixingParametersToDefault();

  // setup integration settings
  nus.Set_h_max( 500.0*nus.units.km );
  nus.Set_rel_error(1.0e-12);
  nus.Set_abs_error(1.0e-12);

  marray<double,1> E_range = nus.GetERange();

  // construct the initial state
  marray<double,2> inistate{60,3};
  double N0 = 1.0e18;
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        // initialze muon state
        inistate[i][k] = N0*pow(E_range[i],-1.0);
      }
  }

  // set the initial state
  nus.Set_initial_state(inistate,flavor);

  nus.EvolveState();

  nus.WriteStateHDF5("./hdf5_muon_test.hdf5");

  nuSQUIDS nus_read("./hdf5_muon_test.hdf5");

  //===================== TESTS ==========================//
  //===================== TESTS ==========================//
  //===================== TESTS ==========================//

  // checking that times are the same
  double squids_time_diff = nus.Get_t() - nus_read.Get_t();
  if (squids_time_diff > 1.0e-15)
    std::cout << nus.Get_t() << " " << nus_read.Get_t() << std::endl;

  double squids_time_initial_diff = nus.Get_t_initial() - nus_read.Get_t_initial();
  if (squids_time_initial_diff > 1.0e-15)
    std::cout << nus.Get_t_initial() << " " << nus_read.Get_t_initial() << std::endl;

  // check that number of neutrinos are the same
  if ( nus.GetNumNeu() != nus_read.GetNumNeu())
    std::cout << nus.GetNumNeu() << " " << nus_read.GetNumNeu() << std::endl;

  // check that the energy ranges are the same
  if ( nus.GetNumE() != nus_read.GetNumE() )
    std::cout << nus.GetNumE() << " " << nus_read.GetNumE() << std::endl;

  // check that the energies are the same
  marray<double,1> energy_diff = nus.GetERange() - nus_read.GetERange();
  for (int i = 0 ; i < nus_read.GetNumE(); i++){
    if ( energy_diff[i] > 1.0e-15)
      std::cout << "ED " << i << " " << energy_diff[i] << std::endl;
  }

  // checking that the mixing angles are the same
  for ( int i = 0 ; i < nus.GetNumNeu() ; i++){
    for ( int j = 0; j < i ; j ++){
      double diff_th = nus.Get_MixingAngle(j,i) - nus_read.Get_MixingAngle(j,i);
      double diff_cp = nus.Get_CPPhase(j,i) - nus_read.Get_CPPhase(j,i);
      if ( diff_th > 1.0e-15 )
        std::cout << "MA " << nus.Get_MixingAngle(j,i) << " " << nus_read.Get_MixingAngle(j,i) << std::endl;
      if ( diff_cp > 1.0e-15 )
        std::cout << "CP " << nus.Get_CPPhase(j,i) << " " << nus_read.Get_CPPhase(j,i) << std::endl;
    }
  }

  // check that square mass differences are the same
  for ( int i = 1; i < nus.GetNumNeu() ; i++){
    double diff_dmsq = nus.Get_SquareMassDifference(i) - nus_read.Get_SquareMassDifference(i);
    if (diff_dmsq > 1.0e-15)
      std::cout << "DMSQ "<< nus.Get_SquareMassDifference(i) << " " << nus_read.Get_SquareMassDifference(i) << std::endl;
  }

  // check that the projectors are the same
  for ( int i = 0 ; i < nus_read.GetNumNeu(); i++){
    squids::SU_vector fproj_diff = nus.GetFlavorProj(i) - nus_read.GetFlavorProj(i);
    for ( double component : fproj_diff.GetComponents()){
      if (component > 1.0e-15)
        std::cout << "FP" << component << std::endl;
    }
    squids::SU_vector mproj_diff = nus.GetMassProj(i) - nus_read.GetMassProj(i);
    for ( double component : mproj_diff.GetComponents()){
      if (component > 1.0e-15)
        std::cout << "MP" << component << std::endl;
    }
  }

  // check that the body is the same
  if ( nus.GetBody()->GetId() != nus_read.GetBody()->GetId() )
    std::cout << "DiffBody " << nus.GetBody()->GetId() << " " << nus_read.GetBody()->GetId() << std::endl;

  // check that the body parameters are the same
  if ( nus.GetBody()->GetBodyParams().size() != nus.GetBody()->GetBodyParams().size())
    std::cout << "BodyParamsSize " << nus.GetBody()->GetBodyParams().size() << " " << nus_read.GetBody()->GetBodyParams().size() << std::endl;

  for ( int i = 0 ; i < nus.GetBody()->GetBodyParams().size(); i++){
    double body_params_diff = nus.GetBody()->GetBodyParams()[i] - nus_read.GetBody()->GetBodyParams()[i];
    if (body_params_diff > 1.0e-15)
      std::cout << "BP " << body_params_diff << std::endl;
  }

  // check that track parameters are the same
  if ( nus.GetTrack()->GetTrackParams().size() != nus.GetTrack()->GetTrackParams().size())
    std::cout << "TrackParamsSize " << nus.GetTrack()->GetTrackParams().size() << " " << nus_read.GetTrack()->GetTrackParams().size() << std::endl;

  for ( int i = 0 ; i < nus.GetTrack()->GetTrackParams().size(); i++){
    double track_params_diff = nus.GetTrack()->GetTrackParams()[i] - nus_read.GetTrack()->GetTrackParams()[i];
    if (track_params_diff > 1.0e-15)
      std::cout << "TP " << track_params_diff << std::endl;
  }

  double diff_x_current = nus.GetTrack()->GetX() - nus_read.GetTrack()->GetX();
  if ( diff_x_current > 1.0e-15 )
    std::cout << "XC " << nus.GetTrack()->GetX() << " " << nus_read.GetTrack()->GetX() << std::endl;

  double diff_x_initial = nus.GetTrack()->GetInitialX() - nus_read.GetTrack()->GetInitialX();
  if ( diff_x_initial > 1.0e-15 )
    std::cout << "XI " << nus.GetTrack()->GetInitialX() << " " << nus_read.GetTrack()->GetInitialX() << std::endl;

  double diff_x_final = nus.GetTrack()->GetFinalX() - nus_read.GetTrack()->GetFinalX();
  if ( diff_x_final > 1.0e-15 )
    std::cout << "XF " << nus.GetTrack()->GetFinalX() << " " << nus_read.GetTrack()->GetFinalX() << std::endl;

  // checking that hamiltonian is the same
  for ( int i = 0 ; i < nus_read.GetNumE(); i++){
    squids::SU_vector hamiltonian_diff = nus.GetHamiltonian(i) - nus_read.GetHamiltonian(i);
    for (double component : hamiltonian_diff.GetComponents()){
      if (component > 1.0e-15)
        std::cout << "HC" << component << std::endl;
    }
  }

  // check that the state is the same
  for ( int i = 0 ; i < nus_read.GetNumE(); i++){
    squids::SU_vector state_diff = nus.GetState(i) - nus_read.GetState(i);
    for (double component : state_diff.GetComponents()){
      if (component > 1.0e-15)
        std::cout << "SC" << component << std::endl;
    }
  }

  // check that the expectation values are the same
  for ( int i = 0 ; i < nus_read.GetNumE(); i++){
      for ( int k = 0; k < nus_read.GetNumNeu(); k++){
        double dif = fabs(nus.EvalFlavorAtNode(k,i) - nus_read.EvalFlavorAtNode(k,i));
        if (dif > 1.0e-15)
          std::cout << "DIF" << i << " " << k << " " << dif << std::endl;
      }
  }

  return 0;
}
