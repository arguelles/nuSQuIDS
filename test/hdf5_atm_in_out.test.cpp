#include <nuSQuIDS/nuSQuIDS.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace nusquids;

int main(){

  squids::Const units;
  nuSQUIDSAtm<nuSQUIDS> nus(-1,0,5,1.e2*units.GeV,1.e6*units.GeV,60,3,neutrino,true,true);

  double phi = acos(-0.5);
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track_atm = std::make_shared<EarthAtm::Track>(phi);

  // set mixing angles and masses
  nus.Set_MixingParametersToDefault();

  // setup integration settings
  nus.Set_h_max( 500.0*units.km );
  nus.Set_rel_error(1.0e-12);
  nus.Set_abs_error(1.0e-12);

  marray<double,1> E_range = nus.GetERange();

  // construct the initial state
  marray<double,3> inistate{5,60,3};
  double N0 = 1.0e18;
  for ( int ci = 0 ; ci < inistate.extent(0); ci++){
    for ( int ei = 0 ; ei < inistate.extent(0); ei++){
        for ( int k = 0; k < inistate.extent(1); k ++){
          // initialze muon state
          inistate[ci][ei][k] = (1+ci)*N0*pow(E_range[ei],-1.0);
        }
    }
  }

  // set the initial state
  nus.Set_initial_state(inistate,flavor);

  nus.EvolveState();

  nus.WriteStateHDF5("./hdf5_atm_muon_test.hdf5");

  nuSQUIDSAtm<nuSQUIDS> nus_read("./hdf5_atm_muon_test.hdf5");

  //===================== TESTS ==========================//
  //===================== TESTS ==========================//
  //===================== TESTS ==========================//

  if(nus.GetnuSQuIDS().size() != nus_read.GetnuSQuIDS().size()){
    std::cout << nus.GetnuSQuIDS().size() << " " << nus_read.GetnuSQuIDS().size() << std::endl;
    return 1;
  }

  unsigned int counter = 0;
  // checking that times are the same
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    double squids_time_diff = nsq.Get_t() - nsq_read.Get_t();
    if (std::abs(squids_time_diff) > 1.0e-15)
      std::cout << nsq.Get_t() << " " << nsq_read.Get_t() << std::endl;
    counter++;
  }

  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    double squids_time_initial_diff = nsq.Get_t_initial() - nsq_read.Get_t_initial();
    if (std::abs(squids_time_initial_diff) > 1.0e-15)
      std::cout << nsq.Get_t_initial() << " " << nsq_read.Get_t_initial() << std::endl;
    counter++;
  }

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
      if ( std::abs(diff_th) > 1.0e-15 )
        std::cout << "MA " << nus.Get_MixingAngle(j,i) << " " << nus_read.Get_MixingAngle(j,i) << std::endl;
      if ( std::abs(diff_cp) > 1.0e-15 )
        std::cout << "CP " << nus.Get_CPPhase(j,i) << " " << nus_read.Get_CPPhase(j,i) << std::endl;
    }
  }

  // check that square mass differences are the same
  for ( int i = 1; i < nus.GetNumNeu() ; i++){
    double diff_dmsq = nus.Get_SquareMassDifference(i) - nus_read.Get_SquareMassDifference(i);
    if (std::abs(diff_dmsq) > 1.0e-15)
      std::cout << "DMSQ "<< nus.Get_SquareMassDifference(i) << " " << nus_read.Get_SquareMassDifference(i) << std::endl;
  }

  // check that the projectors are the same
  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    for ( int i = 0 ; i < nus_read.GetNumNeu(); i++){
      squids::SU_vector fproj_diff = nsq.GetFlavorProj(i) - nsq_read.GetFlavorProj(i);
      for ( double component : fproj_diff.GetComponents()){
        if (std::abs(component) > 1.0e-15)
          std::cout << "FP" << component << std::endl;
      }
      squids::SU_vector mproj_diff = nsq.GetMassProj(i) - nsq_read.GetMassProj(i);
      for ( double component : mproj_diff.GetComponents()){
        if (std::abs(component) > 1.0e-15)
          std::cout << "MP" << component << std::endl;
      }
    }
    counter++;
  }

  // check that the body is the same
  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    if ( nsq.GetBody()->GetId() != nsq_read.GetBody()->GetId() )
      std::cout << "DiffBody " << nsq.GetBody()->GetId() << " " << nsq_read.GetBody()->GetId() << std::endl;
    counter++;
  }

  // check that the body parameters are the same
  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    if ( nsq.GetBody()->GetBodyParams().size() != nsq_read.GetBody()->GetBodyParams().size())
      std::cout << "BodyParamsSize " << nsq.GetBody()->GetBodyParams().size() << " " << nsq_read.GetBody()->GetBodyParams().size() << std::endl;
    counter++;
  }

  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    for ( int i = 0 ; i < nsq.GetBody()->GetBodyParams().size(); i++){
      double body_params_diff = nsq.GetBody()->GetBodyParams()[i] - nsq_read.GetBody()->GetBodyParams()[i];
      if (std::abs(body_params_diff) > 1.0e-15)
        std::cout << "BP " << body_params_diff << std::endl;
    }
    counter++;
  }

  // check that track parameters are the same
  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    if ( nsq.GetTrack()->GetTrackParams().size() != nsq_read.GetTrack()->GetTrackParams().size())
      std::cout << "TrackParamsSize " << nsq.GetTrack()->GetTrackParams().size() << " " << nsq_read.GetTrack()->GetTrackParams().size() << std::endl;
    counter++;
  }

  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    for ( int i = 0 ; i < nsq.GetTrack()->GetTrackParams().size(); i++){
      double track_params_diff = nsq.GetTrack()->GetTrackParams()[i] - nsq_read.GetTrack()->GetTrackParams()[i];
      if (std::abs(track_params_diff) > 1.0e-15)
        std::cout << "TP " << track_params_diff << std::endl;
    }
    counter++;
  }

  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    double diff_x_current = nsq.GetTrack()->GetX() - nsq_read.GetTrack()->GetX();
    if ( std::abs(diff_x_current) > 1.0e-15 )
      std::cout << "XC " << nsq.GetTrack()->GetX() << " " << nsq_read.GetTrack()->GetX() << std::endl;
    counter++;
  }

  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    double diff_x_initial = nsq.GetTrack()->GetInitialX() - nsq_read.GetTrack()->GetInitialX();
    if ( std::abs(diff_x_initial) > 1.0e-15 )
      std::cout << "XI " << nsq.GetTrack()->GetInitialX() << " " << nsq_read.GetTrack()->GetInitialX() << std::endl;
    counter++;
  }

  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    double diff_x_final = nsq.GetTrack()->GetFinalX() - nsq_read.GetTrack()->GetFinalX();
    if ( std::abs(diff_x_final) > 1.0e-15 )
      std::cout << "XF " << nsq.GetTrack()->GetFinalX() << " " << nsq_read.GetTrack()->GetFinalX() << std::endl;
    counter++;
  }

  // checking that hamiltonian is the same
  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    for ( int i = 0 ; i < nus_read.GetNumE(); i++){
      squids::SU_vector hamiltonian_diff = nsq.GetHamiltonian(i) - nsq_read.GetHamiltonian(i);
      for (double component : hamiltonian_diff.GetComponents()){
        if (std::abs(component) > 1.0e-15)
          std::cout << "HC" << component << std::endl;
      }
    }
    counter++;
  }

  // check that the state is the same
  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    for ( int i = 0 ; i < nus_read.GetNumE(); i++){
      squids::SU_vector state_diff = nsq.GetState(i) - nsq_read.GetState(i);
      for (double component : state_diff.GetComponents()){
        if (std::abs(component) > 1.0e-15)
          std::cout << "SC" << component << std::endl;
      }
    }
    counter++;
  }

  // check that the expectation values are the same
  counter = 0;
  for (nuSQUIDS &nsq : nus.GetnuSQuIDS()){
    nuSQUIDS &nsq_read = nus_read.GetnuSQuIDS(counter);
    for ( int i = 0 ; i < nsq_read.GetNumE(); i++){
        for ( int k = 0; k < nsq_read.GetNumNeu(); k++){
          double dif = fabs(nsq.EvalFlavorAtNode(k,i) - nsq_read.EvalFlavorAtNode(k,i));
          if (std::abs(dif) > 1.0e-15)
            std::cout << "DIF" << i << " " << k << " " << dif << std::endl;
        }
    }
    counter++;
  }

  marray<double,1> costh_range = nus.GetCosthRange();
  marray<double,1> enu_range = nus.GetERange();

  for ( int ci = 0 ; ci < inistate.extent(0); ci++){
    for ( int ei = 0 ; ei < inistate.extent(1)-1; ei++){
      for ( int k = 0; k < inistate.extent(2); k ++){
        double costh = costh_range[ci];
        double enu = enu_range[ei]+0.1;

        double dif = fabs(nus.EvalFlavor(k,costh,enu) - nus_read.EvalFlavor(k,costh,enu));
        if (std::abs(dif) > 1.0e-15)
          std::cout << "DIF" << costh << " " << enu << " " << k << " " << dif << std::endl;
      }
    }
  }
  return 0;
}
