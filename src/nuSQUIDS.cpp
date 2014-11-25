#include "nuSQUIDS.h"

namespace nusquids{

void nuSQUIDS::init(void){
  // single energy implementation
  ne = 1;

  if (NT == "neutrino" || NT == "antineutrino")
    nrhos = 1;
  else {
    throw std::runtime_error("nuSQUIDS::Error::NT = {neutrino,antineutrino} not : " + NT);
  }

  if ( numneu > 6 )
    throw std::runtime_error("nuSQUIDS::Error::Maximum number of neutrinos exceded");
  //assert( "nuSQUIDS::Error::Maximum number of neutrinos exceded" && numneu < 6);
  nsun = numneu;

  //initialize SQUIDS
  ini(ne,numneu,1,0,0);
  Set_CoherentInteractions(true);
  Set_h_max(std::numeric_limits<double>::max() );

  //===============================
  // set parameters to default   //
  //===============================

  Set_MixingParametersToDefault();

  //===============================
  // physics CP sign for aneu    //
  //===============================
  if ( NT == "antineutrino" ){
    Set(DELTA1,params.GetPhase(param_label_index[DELTA1][0],param_label_index[DELTA1][1]));
    Set(DELTA2,params.GetPhase(param_label_index[DELTA2][0],param_label_index[DELTA2][1]));
    Set(DELTA3,params.GetPhase(param_label_index[DELTA3][0],param_label_index[DELTA3][1]));
  }


  //===============================
  // init projectors             //
  //===============================

  iniProjectors();

  //===============================
  // init square mass difference //
  //===============================

  H0_array.resize(ne);
  for(int ie = 0; ie < ne; ie++){
    H0_array[ie] = SU_vector(nsun);
  }

  iniH0();

  //===============================
  // END                         //
  //===============================
  inusquids = true;
}

void nuSQUIDS::Set_E(double Enu){
  if( ne != 1 )
    throw std::runtime_error("nuSQUIDS::Error:Cannot use Set_E in single energy mode.");
  E_range = std::vector<double>{Enu};
  Set_xrange(E_range[0],E_range[ne-1],"lin");

  ienergy = true;
}

void nuSQUIDS::init(double Emin,double Emax,int Esize)
{

  if (NT == "neutrino" || NT == "antineutrino")
    nrhos = 1;
  else if (NT == "both")
    nrhos = 2;
  else {
    throw std::runtime_error("nuSQUIDS::Error::NT = {neutrino,antineutrino,both} not : " + NT);
  }

  if ( numneu > 6 )
    throw std::runtime_error("nuSQUIDS::Error::Maximum number of neutrinos exceded");
  nsun = numneu;
  if ( Emax < Emin )
    throw std::runtime_error("nuSQUIDS::Error::Emax < Emin.");

  ne = Esize;
  if ( Esize <= 0 )
    throw std::runtime_error("nuSQUIDS::Error::Esize must be at greater than zero.");

  //===============================
  // BEGIN                       //
  //===============================

  // initialize SQUIDS
  if (iinteraction)
    ini(ne,numneu,nrhos,nrhos,0);
  else
    ini(ne,numneu,nrhos,0,0);

  SetScalarsToZero();

  t = 0;

  Set_CoherentInteractions(true);
  Set_h_max(std::numeric_limits<double>::max());

  //===============================
  // initialize energy arrays    //
  //===============================
  if(elogscale){
    E_range = logspace(Emin*params.GeV,Emax*params.GeV,ne-1);
    Set_xrange(E_range[0],E_range[Esize-1],"log");
  }
  else{
    E_range = linspace(Emin*params.GeV,Emax*params.GeV,ne-1);
    Set_xrange(E_range[0],E_range[Esize-1],"lin");
  }
  delE.resize(ne-1);
  for(int ei = 0; ei < ne -1; ei++)
    delE[ei] = E_range[ei+1] - E_range[ei];

  ienergy = true;

  //===============================
  // set parameters to default   //
  //===============================

  Set_MixingParametersToDefault();

  //===============================
  // init projectors             //
  //===============================

  iniProjectors();

  //===============================
  // init square mass difference //
  //===============================

  H0_array.resize(ne);
  for(int ie = 0; ie < ne; ie++){
    H0_array[ie] = SU_vector(nsun);
  }

  iniH0();

  if(iinteraction){
    //===============================
    // init XS and TDecay objects  //
    //===============================

    // initialize cross section object
    ncs.Init(E_range[0],E_range[ne-1],ne-1);
    // initialize tau decay spectra object
    tdc.Init(E_range[0],E_range[ne-1],ne-1);
    // initialize cross section and interaction arrays
    dNdE_NC.resize(nrhos); dNdE_CC.resize(nrhos);
    invlen_NC.resize(nrhos); invlen_CC.resize(nrhos); invlen_INT.resize(nrhos);
    sigma_CC.resize(nrhos); sigma_NC.resize(nrhos);
    for(int rho = 0; rho < nrhos; rho++){
      dNdE_NC[rho].resize(numneu); dNdE_CC[rho].resize(numneu);
      invlen_NC[rho].resize(numneu); invlen_CC[rho].resize(numneu); invlen_INT[rho].resize(numneu);
      sigma_CC[rho].resize(numneu); sigma_NC[rho].resize(numneu);
      for(int flv = 0; flv < numneu; flv++){
          dNdE_NC[rho][flv].resize(ne); dNdE_CC[rho][flv].resize(ne);
          invlen_NC[rho][flv].resize(ne); invlen_CC[rho][flv].resize(ne); invlen_INT[rho][flv].resize(ne);
          sigma_CC[rho][flv].resize(ne); sigma_NC[rho][flv].resize(ne);
          for(int e1 = 0; e1 < ne; e1++){
            dNdE_NC[rho][flv][e1].resize(e1); dNdE_CC[rho][flv][e1].resize(e1);
          }
      }
    }
    // initialize the tau decay and interaction array
    invlen_tau.resize(ne);
    dNdE_tau_all.resize(ne); dNdE_tau_lep.resize(ne);
    for(int e1 = 0; e1 < ne; e1++){
      dNdE_tau_all[e1].resize(e1); dNdE_tau_lep[e1].resize(e1);
    }
  }
  //===============================
  // Tau properties              //
  //===============================

  taubr_lep = 0.14;
  tau_lifetime = 2.906e-13*params.sec;
  tau_mass = 1776.82*params.MeV;
  tau_reg_scale = 100.0*params.km;

  //===============================
  // Fill in arrays              //
  //===============================
  if(iinteraction){
    InitializeInteractions();
    Set_NonCoherentInteractions(true);
    Set_ScalarInteractions(true);
    Set_OtherInteractions(true);
  }

  //===============================
  // END                         //
  //===============================

  inusquids = true;
}

void nuSQUIDS::PreDerive(double x){
  track->SetX(tunit*x);
  if( basis != mass){
    EvolveProjectors(tunit*x);
  }
  if(iinteraction){
    UpdateInteractions();
  }
  if(progressbar and progressbar_count%progressbar_loop ==0 ){
    ProgressBar();
  }
  progressbar_count++;
  AddToPreDerive(x);
}

void nuSQUIDS::EvolveProjectors(double x){
  for(int rho = 0; rho < nrhos; rho++){
    for(int flv = 0; flv < numneu; flv++){
      for(int ei = 0; ei < ne; ei++){
        // will only evolve the flavor projectors
        //evol_b0_proj[rho][flv][ei] = b0_proj[flv].SUEvolve(h0,(x-t_ini));
        evol_b1_proj[rho][flv][ei] = b1_proj[rho][flv].SUEvolve(H0_array[ei],(x-t_ini));
      }
    }
  }
}

SU_vector nuSQUIDS::H0(double Enu){
  return DM2*(0.5/Enu);
}

SU_vector nuSQUIDS::HI(int ei){
    double ye = body->ye(track);
    double density = body->density(track);

    double CC = params.sqrt2*params.GF*params.Na*pow(params.cm,-3)*density*ye;
    double NC;

    if (ye < 1.0e-10){
      NC = params.sqrt2*params.GF*params.Na*pow(params.cm,-3)*density;
    }
    else {
      NC = CC*(-0.5*(1.0-ye)/ye);
    }

    // construct potential in flavor basis
    SU_vector potential = (CC+NC)*evol_b1_proj[index_rho][0][ei];
    potential += (NC)*(evol_b1_proj[index_rho][1][ei]);
    potential += (NC)*(evol_b1_proj[index_rho][2][ei]);

    if (basis == mass){
      potential += H0_array[ei];
    }

    if ((index_rho == 0 and NT=="both") or NT=="neutrino"){
        // neutrino potential
        return potential;
    } else if ((index_rho == 1 and NT=="both") or NT=="antineutrino"){
        // antineutrino potential
        return (-1.0)*potential;
    } else{
        throw std::runtime_error("nuSQUIDS::HI : unknown particle or antiparticle");
    }
}

SU_vector nuSQUIDS::GammaRho(int ei){
    SU_vector V(nsun);
    if (not iinteraction){
      return V;
    }

    //std::cout << invlen_INT[index_rho][0][ei] << " " << invlen_INT[index_rho][0][ei] << " " << invlen_INT[index_rho][0][ei] << std::endl;
    //std::cout << "GammaRho" << std::endl;
    //std::cout << evol_b1_proj[index_rho][0][ei]*(0.5*invlen_INT[index_rho][0][ei]) +
    //       evol_b1_proj[index_rho][1][ei]*(0.5*invlen_INT[index_rho][1][ei]) +
    //       evol_b1_proj[index_rho][2][ei]*(0.5*invlen_INT[index_rho][2][ei]) << std::endl;


    V = evol_b1_proj[index_rho][0][ei]*(0.5*invlen_INT[index_rho][0][ei]);
    V += evol_b1_proj[index_rho][1][ei]*(0.5*invlen_INT[index_rho][1][ei]);
    V += evol_b1_proj[index_rho][2][ei]*(0.5*invlen_INT[index_rho][2][ei]);

    return V;
}

SU_vector nuSQUIDS::InteractionsRho(int e1){
  SU_vector nc_term(nsun);

  if (not iinteraction){
    return nc_term;
  }

  // this implements the NC interactinos
  // the tau regeneration terms are implemented at the end
  SU_vector temp1, temp2;
  for(int e2 = e1 + 1; e2 < ne; e2++){
    // here we assume the cross section to be the same for all flavors
    //std::cout << dNdE_NC[index_rho][0][e2][e1] << " " << invlen_NC[index_rho][0][e2] << std::endl;
    temp1 = evol_b1_proj[index_rho][0][e1] + evol_b1_proj[index_rho][1][e1];
    temp1 += evol_b1_proj[index_rho][2][e1];
    temp2 = ACommutator(temp1,state[e2].rho[index_rho]);
    nc_term += temp2*(0.5*dNdE_NC[index_rho][0][e2][e1]*invlen_NC[index_rho][0][e2]);
  }

  //std::cout << "nc term" << std::endl;
  //std::cout << nc_term << std::endl;

  return nc_term;
}

double nuSQUIDS::GammaScalar(int ei){
  // we will just keep all the taus and convert them at the end
  return 0.0;
}

double nuSQUIDS::InteractionsScalar(int ei){
  if (not iinteraction){
    return 0.0;
  }

  double nutautoleptau = 0.0;
  for(int e2 = ei + 1; e2 < ne; e2++)
    nutautoleptau += (evol_b1_proj[index_scalar][2][ei]*state[e2].rho[index_scalar])*
                     (invlen_CC[index_scalar][2][ei])*(dNdE_CC[index_scalar][2][e2][ei])*delE[e2];
  return nutautoleptau;
}

double nuSQUIDS::GetNucleonNumber(){
    double density = body->density(track);
    double num_nuc = (params.gr*pow(params.cm,-3))*density*2.0/(params.proton_mass+params.neutron_mass);

    #ifdef UpdateInteractions_DEBUG
        cout << "Density " << density << endl;
        cout << "Nucleon Number " << num_nuc << endl;
    #endif

    if(num_nuc == 0){
        num_nuc = params.Na*pow(params.cm,-3)*1.0e-10;
    }

    return num_nuc;
}

void nuSQUIDS::UpdateInteractions(){
    double num_nuc = GetNucleonNumber();
    for(int rho = 0; rho < nrhos; rho++){
      for(int flv = 0; flv < numneu; flv++){
              #ifdef UpdateInteractions_DEBUG
              cout << "============" << flv << "============" << endl;
              #endif
          for(int e1 = 0; e1 < ne; e1++){
              #ifdef UpdateInteractions_DEBUG
                  cout << "== CC NC Terms x = " << track->x/params.km << " [km] ";
                  cout << "E = " << x[e1] << " [eV] ==" << endl;
                  cout << "CC : " << sigma_CC[rho][flv][e1]*num_nuc << " NC : " << sigma_NC[rho][flv][e1]*num_nuc << endl;
                  cout << "==" << endl;
              #endif
              invlen_NC[rho][flv][e1] = sigma_NC[rho][flv][e1]*num_nuc;
              invlen_CC[rho][flv][e1] = sigma_CC[rho][flv][e1]*num_nuc;
              invlen_INT[rho][flv][e1] = invlen_NC[rho][flv][e1] + invlen_CC[rho][flv][e1];
          }
      }
    }
}

void nuSQUIDS::InitializeInteractions(){

    //units
    double cm2GeV = pow(params.cm,2)*pow(params.GeV,-1);
    double cm2 = pow(params.cm,2);
    double GeVm1 = pow(params.GeV,-1);

    // load cross sections
    // initializing cross section arrays temporary array
    array4D dsignudE_CC,dsignudE_NC;

    // filling cross section arrays
    dsignudE_NC.resize(nrhos); dsignudE_CC.resize(nrhos);
    for(int neutype = 0; neutype < nrhos; neutype++){
      dsignudE_CC[neutype].resize(numneu); dsignudE_NC[neutype].resize(numneu);
      for(int flv = 0; flv < numneu; flv++){
          dsignudE_CC[neutype][flv].resize(ne); dsignudE_NC[neutype][flv].resize(ne);
          for(int e1 = 0; e1 < ne; e1++){
              dsignudE_CC[neutype][flv][e1].resize(e1); dsignudE_NC[neutype][flv][e1].resize(e1);
              // differential cross sections
              for(int e2 = 0; e2 < e1; e2++){
                  dsignudE_NC[neutype][flv][e1][e2] = ncs.dsde_NC(e1,e2,flv,neutype)*cm2GeV;
                  dsignudE_CC[neutype][flv][e1][e2] = ncs.dsde_CC(e1,e2,flv,neutype)*cm2GeV;
              }
              // total cross sections
              sigma_CC[neutype][flv][e1] = ncs.sigma_CC(e1,flv,neutype)*cm2;
              sigma_NC[neutype][flv][e1] = ncs.sigma_NC(e1,flv,neutype)*cm2;
          }
      }
    }

    #ifdef FixCrossSections
    // fix charge current and neutral current differential cross sections
    for(int neutype = 0; neutype < nrhos; neutype++){
      double XCC_MIN,XNC_MIN,XCC_int,XNC_int,CC_rescale,NC_rescale;
      for(int flv = 0; flv < numneu; flv++){
          XCC_MIN = sigma_CC[neutype][flv][0];
          XNC_MIN = sigma_CC[neutype][flv][0];
          for(int e1 = 0; e1 < ne; e1++){
              XCC_int = 0.0;
              XNC_int = 0.0;
              for(int e2 = 0; e2 < e1; e2++){
                  XCC_int += dsignudE_CC[neutype][flv][e1][e2]*delE[e2];
                  XNC_int += dsignudE_NC[neutype][flv][e1][e2]*delE[e2];
              }

              if(e1 != 0 ){
                  CC_rescale = (sigma_CC[neutype][flv][e1] - XCC_MIN)/XCC_int;
                  NC_rescale = (sigma_NC[neutype][flv][e1] - XNC_MIN)/XNC_int;

                  for(int e2 = 0; e2 < e1; e2++){
                      dsignudE_CC[neutype][flv][e1][e2] = dsignudE_CC[neutype][flv][e1][e2]*CC_rescale;
                      dsignudE_NC[neutype][flv][e1][e2] = dsignudE_NC[neutype][flv][e1][e2]*NC_rescale;
                  }
              }
          }
      }
    }
    #endif

    // constructing dNdE
    for(int rho = 0; rho < nrhos; rho++){
      for(int flv = 0; flv < numneu; flv++){
          for(int e1 = 0; e1 < ne; e1++){
              for(int e2 = 0; e2 < e1; e2++){
                  if (dsignudE_NC[rho][flv][e1][e2] < 1.0e-50 or (dsignudE_NC[rho][flv][e1][e2] != dsignudE_NC[rho][flv][e1][e2])){
                      dNdE_NC[rho][flv][e1][e2] = 0.0;
                  } else {
                      dNdE_NC[rho][flv][e1][e2] = (dsignudE_NC[rho][flv][e1][e2])/(sigma_NC[rho][flv][e1]);
                  }
                  if (dsignudE_CC[rho][flv][e1][e2] < 1.0e-50 or (dsignudE_CC[rho][flv][e1][e2] != dsignudE_CC[rho][flv][e1][e2])){
                      dNdE_CC[rho][flv][e1][e2] = 0.0;
                  } else {
                      dNdE_CC[rho][flv][e1][e2] = (dsignudE_CC[rho][flv][e1][e2])/(sigma_CC[rho][flv][e1]);
                  }
              }
          }
      }
    }

    // initialize interaction lenghts to zero
    // tau decay length array
    for(int e1 = 0; e1 < ne; e1++){
        invlen_tau[e1] = 1.0/(tau_lifetime*E_range[e1]*tau_mass);
    }

    // load tau decay spectra

    // constructing dNdE_tau_lep/dNdE_tau_all
    for(int e1 = 0; e1 < ne; e1++){
        for(int e2 = 0; e2 < e1; e2++){
            dNdE_tau_all[e1][e2] = tdc.dNdEnu_All(e1,e2)*GeVm1;
            dNdE_tau_lep[e1][e2] = tdc.dNdEnu_Lep(e1,e2)*GeVm1;
        }
    }

    #ifdef FixCrossSections
    // fix tau decay spectra cross section
    double tau_all_int,tau_lep_int,tau_lep_rescale,tau_all_rescale;
    for(int e1 = 1; e1 < ne; e1++){
        tau_all_int = 0.0;
        tau_lep_int = 0.0;
        for(int e2 = 0; e2 < e1; e2++){
             tau_all_int += dNdE_tau_all[e1][e2]*delE[e2];
             tau_lep_int += dNdE_tau_lep[e1][e2]*delE[e2];
        }

        if( dNdE_tau_all[e1][0]*E_range[0] < 0.25 ) {
            tau_all_rescale = (1.0 - dNdE_tau_all[e1][0]*E_range[0])/tau_all_int;
            tau_lep_rescale = (taubr_lep - dNdE_tau_lep[e1][0]*E_range[0])/tau_lep_int;

            for(int e2 = 0; e2 < e1; e2++){
                dNdE_tau_all[e1][e2] = dNdE_tau_all[e1][e2]*tau_all_rescale;
                dNdE_tau_lep[e1][e2] = dNdE_tau_lep[e1][e2]*tau_lep_rescale;
            }
        }
    }
    #endif
}

void nuSQUIDS::Set_Body(std::shared_ptr<Body> body_in){
  body = body_in;
  ibody = true;
}

void nuSQUIDS::Set_Track(std::shared_ptr<Track> track_in){
  track = track_in;
  itrack = true;
}

void nuSQUIDS::EvolveState(){
  // check for BODY and TRACK status
  if ( body == NULL )
    throw std::runtime_error("nuSQUIDS::Error::BODY is a NULL pointer");
  if (not ibody )
    throw std::runtime_error("nuSQUIDS::Error::BODY is a NULL pointer");
  if ( track == NULL )
    throw std::runtime_error("nuSQUIDS::Error::TRACK is a NULL pointer");
  if ( not itrack )
    throw std::runtime_error("nuSQUIDS::Error::TRACK is not initialized");
  if ( not istate )
    throw std::runtime_error("nuSQUIDS::Error::Initial state not initialized");
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");

  // restart the clock
  Set_Initial_Time();
  // initializing the projectors and hamiltonian
  SetIniFlavorProyectors();
  iniH0();

  if( not tauregeneration ){
    EvolveSUN(track->GetFinalX()-track->GetInitialX());
  }
  else {
    int tau_steps = (int) ((track->GetFinalX() - track->GetInitialX())/tau_reg_scale);
    std::cout << tau_steps << " " << tau_reg_scale/units.km<< std::endl;
    for (int i = 0; i < tau_steps; i++){
      double x_inter = track->GetInitialX() + (double)(i)*tau_reg_scale;
      std::cout << x_inter/units.km << std::endl;
      EvolveSUN(tau_reg_scale);
      ConvertTauIntoNuTau();
    }
    EvolveSUN(track->GetFinalX()-tau_reg_scale*tau_steps);
    ConvertTauIntoNuTau();
  }
}

void nuSQUIDS::SetScalarsToZero(void){
  for(int rho = 0; rho < nscalars; rho++){
    for(int e1 = 0; e1 < ne; e1++){
        state[e1].scalar[rho] = 0.0;
    }
  }
}

void nuSQUIDS::ConvertTauIntoNuTau(void){

  for(int e1 = 0; e1 < ne; e1++){
      double tau_neu_all  = 0.0;
      double tau_neu_lep  = 0.0;
      double tau_aneu_all = 0.0;
      double tau_aneu_lep = 0.0;

      for(int e2 = e1 +1; e2 < ne; e2++){
          //std::cout << dNdE_tau_all[e2][e1] << " " << delE[e2] << " " << state[e2].scalar[0] << std::endl;
          tau_neu_all  += dNdE_tau_all[e2][e1]*delE[e2]*state[e2].scalar[0];
          tau_neu_lep  += dNdE_tau_lep[e2][e1]*delE[e2]*state[e2].scalar[0];
          tau_aneu_all += dNdE_tau_all[e2][e1]*delE[e2]*state[e2].scalar[1];
          tau_aneu_lep += dNdE_tau_lep[e2][e1]*delE[e2]*state[e2].scalar[1];
      }
      // note that the br_lepton is already included in dNdE_tau_lep
      // adding new fluxes
      //std::cout << tau_neu_all << " " << tau_aneu_all << " " << tau_neu_lep << " " << tau_aneu_lep << std::endl;
      state[e1].rho[0]  += tau_neu_all*evol_b1_proj[0][2][e1] +
                            tau_aneu_lep*evol_b1_proj[0][0][e1] +
                            tau_aneu_lep*evol_b1_proj[0][1][e1];
      state[e1].rho[1]  += tau_aneu_all*evol_b1_proj[1][2][e1] +
                            tau_neu_lep*evol_b1_proj[1][0][e1] +
                            tau_neu_lep*evol_b1_proj[1][1][e1];
  }

  // clean all lepton arrays
  SetScalarsToZero();

}

void nuSQUIDS::Set_Initial_Time(){
  t_ini = 0.0;
  t = 0.0;
}

void nuSQUIDS::Set_initial_state(std::vector<double> v, std::string basis){
  if( v.size() == 0 )
    throw std::runtime_error("nuSQUIDS::Error:Null size input array.");
  if( v.size() != numneu )
    throw std::runtime_error("nuSQUIDS::Error::Initial state size not compatible with number of flavors");
  if( not (basis == "flavor" || basis == "mass" ))
    throw std::runtime_error("nuSQUIDS::Error::BASIS can be : flavor or mass.");
  if( NT == "both" )
    throw std::runtime_error("nuSQUIDS::Error::Only supplied neutrino/antineutrino initial state, but set to both.");
  if( ne != 1 )
    throw std::runtime_error("nuSQUIDS::Error::nuSQUIDS initialized in multienergy mode, while state is only single energy");

  Set_Initial_Time();
  SetIniFlavorProyectors();

  for(int i = 0; i < ne; i++){
    for(int r = 0; r < nrhos; r++){
      if (basis == "flavor"){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(int j = 0; j < v.size(); j++)
        {
          state[i].rho[r] += v[j]*b1_proj[r][j];
        }
      }
      else if (basis == "mass"){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(int j = 0; j < v.size(); j++)
        {
          state[i].rho[r] += v[j]*b0_proj[j];
        }
      }
    }
  }

  istate = true;
};

void nuSQUIDS::Set_initial_state(array2D v, std::string basis){
  if( v.size() == 0 )
    throw std::runtime_error("nuSQUIDS::Error:Null size input array.");
  if( v.size() != ne )
    throw std::runtime_error("nuSQUIDS::Error:Input vector with wrong dimensions.");
  if( not (basis == "flavor" || basis == "mass" ))
    throw std::runtime_error("nuSQUIDS::Error::BASIS can be : flavor or mass.");
  if( NT == "both" )
    throw std::runtime_error("nuSQUIDS::Error::Only supplied neutrino/antineutrino initial state, but set to both.");

  Set_Initial_Time();
  SetIniFlavorProyectors();

  for(int i = 0; i < ne; i++){
    for(int r = 0; r < nrhos; r++){
      if (basis == "flavor"){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(int j = 0; j < numneu; j++)
        {
          state[i].rho[r] += v[i][j]*b1_proj[r][j];
        }
      }
      else if (basis == "mass"){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(int j = 0; j < numneu; j++){
          state[i].rho[r] += v[i][j]*b0_proj[j];
        }
      }
    }
  }

  istate = true;
}

void nuSQUIDS::Set_initial_state(array3D v, std::string basis){
  if( v.size() == 0 )
    throw std::runtime_error("nuSQUIDS::Error:Null size input array.");
  if( v.size() != ne )
    throw std::runtime_error("nuSQUIDS::Error:Input vector with wrong dimensions.");
  if( not (basis == "flavor" || basis == "mass" ))
    throw std::runtime_error("nuSQUIDS::Error::BASIS can be : flavor or mass.");
  if( NT != "both" )
    throw std::runtime_error("nuSQUIDS::Error::Supplied neutrino and antineutrino initial state, but not set to both.");

  Set_Initial_Time();
  SetIniFlavorProyectors();

  for(int i = 0; i < ne; i++){
    for(int r = 0; r < nrhos; r++){
      if (basis == "flavor"){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(int j = 0; j < numneu; j++)
        {
          state[i].rho[r] += v[i][r][j]*b1_proj[r][j];
        }
      }
      else if (basis == "mass"){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(int j = 0; j < numneu; j++){
          state[i].rho[r] += v[i][r][j]*b0_proj[j];
        }
      }
    }
  }
  istate = true;
}

void nuSQUIDS::WriteState(std::string filename){
  //assert("nuSQUIDS::Error::State not initialized." && state != NULL);
  //assert("nuSQUIDS::Error::nuSQUIDS not initialized." && inusquids);
  if(state == NULL)
    throw std::runtime_error("nuSQUIDS::Error::State not initialized.");
  if(not inusquids)
    throw std::runtime_error("nuSQUIDS::Error::nuSQUIDS not initialized.");

  array2D tbl_state;
  for(int ie = 0; ie < ne; ie++){
    for(int rho = 0; rho < nrhos; rho++){
      array1D row;
      row.push_back(E_range[ie]/units.GeV);
      for(int i = 0; i < numneu*numneu; i ++)
        row.push_back(state[ie].rho[rho][i]);
      tbl_state.push_back(row);
    }
  }
  quickwrite(filename,tbl_state);
}

void nuSQUIDS::ReadState(std::string filename){

}

array1D nuSQUIDS::GetERange(void){
  return E_range;
}

size_t nuSQUIDS::GetNumE(void){
  return ne;
}

double nuSQUIDS::EvalMass(int flv,double EE, int rho){
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != "both" )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  if ( basis == mass )
    throw std::runtime_error("nuSQUIDS::Error::Use EvalMassAtNode. Interpolation is not recommended on this basis.");
  return GetExpectationValueD(b0_proj[flv], rho, EE);
}

double nuSQUIDS::EvalFlavor(int flv,double EE, int rho){
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != "both" )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  if ( basis == mass )
    throw std::runtime_error("nuSQUIDS::Error::Use EvalMassAtNode. Interpolation is not recommended on this basis.");
  return GetExpectationValueD(b1_proj[rho][flv], rho, EE);
}

double nuSQUIDS::EvalMassAtNode(int flv, int ei, int rho){
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != "both" )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  if(basis == mass)
    return b0_proj[flv]*state[ei].rho[rho];
  return GetExpectationValue(b0_proj[flv], rho, ei);
}

double nuSQUIDS::EvalFlavorAtNode(int flv, int ei, int rho){
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != "both" )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  if(basis == mass)
    return b1_proj[rho][flv]*state[ei].rho[rho];
  return GetExpectationValue(b1_proj[rho][flv], rho, ei);
}

double nuSQUIDS::EvalMass(int flv){
  if(state == NULL)
    throw std::runtime_error("nuSQUIDS::Error::State not initialized.");
  if(not inusquids)
    throw std::runtime_error("nuSQUIDS::Error::nuSQUIDS not initialized.");
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if( ne != 1 )
    throw std::runtime_error("nuSQUIDS::Error::Use this function only in single energy mode.");
  if( flv >= nsun )
    throw std::runtime_error("nuSQUIDS::Error::Flavor index greater than number of initialized flavors.");

  if(basis == mass)
    return b0_proj[flv]*state[0].rho[0];
  return GetExpectationValue(b0_proj[flv], 0, 0);
}

double nuSQUIDS::EvalFlavor(int flv){
  if(state == NULL)
    throw std::runtime_error("nuSQUIDS::Error::State not initialized.");
  if(not inusquids)
    throw std::runtime_error("nuSQUIDS::Error::nuSQUIDS not initialized.");
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if( ne != 1 )
    throw std::runtime_error("nuSQUIDS::Error::Use this function only in single energy mode.");
  if( flv >= nsun )
    throw std::runtime_error("nuSQUIDS::Error::Flavor index greater than number of initialized flavors.");

  if(basis == mass)
    return b1_proj[0][flv]*state[0].rho[0];
  return GetExpectationValue(b1_proj[0][flv], 0, 0);
}

void nuSQUIDS::iniH0(){
  DM2 = SU_vector(nsun);
  for(int i = 1; i < nsun; i++){
      DM2 += (b0_proj[i])*params.GetSquaredEnergyDifference(i);
  }

  if(ienergy){
    for(int ei = 0; ei < ne; ei++){
      H0_array[ei] = H0(E_range[ei]);
    }
  }
}

void nuSQUIDS::AntineutrinoCPFix(int rho){
    if(NT == "antineutrino" or (NT == "both" and rho == 1)){
      Set(DELTA1,-params.GetPhase(0,1));
      Set(DELTA2,-params.GetPhase(0,2));
      Set(DELTA3,-params.GetPhase(0,3));
    }
}

void nuSQUIDS::iniProjectors(){

  b0_proj.resize(numneu);
  for(int flv = 0; flv < numneu; flv++){
    b0_proj[flv] = SU_vector::Projector(nsun,flv);
  }

  b1_proj.resize(nrhos);
  for(int rho = 0; rho < nrhos; rho++){
    b1_proj[rho].resize(numneu);
    for(int flv = 0; flv < numneu; flv++){
      b1_proj[rho][flv] = SU_vector::Projector(nsun,flv);

      AntineutrinoCPFix(rho);
      b1_proj[rho][flv].RotateToB1(params);
      AntineutrinoCPFix(rho);
    }
  }

  evol_b0_proj.resize(nrhos);
  evol_b1_proj.resize(nrhos);
  for(int rho = 0; rho < nrhos; rho++){
    evol_b0_proj[rho].resize(numneu);
    evol_b1_proj[rho].resize(numneu);
    for(int flv = 0; flv < numneu; flv++){
      evol_b0_proj[rho][flv].resize(ne);
      evol_b1_proj[rho][flv].resize(ne);
      for(int e1 = 0; e1 < ne; e1++){
        evol_b0_proj[rho][flv][e1] = SU_vector::Projector(nsun,flv);
        evol_b1_proj[rho][flv][e1] = SU_vector::Projector(nsun,flv);

        AntineutrinoCPFix(rho);
        evol_b1_proj[rho][flv][e1].RotateToB1(params);
        AntineutrinoCPFix(rho);
      }
    }
  }

}

void nuSQUIDS::SetIniFlavorProyectors(){
  for(int rho = 0; rho < nrhos; rho++){
    for(int flv = 0; flv < numneu; flv++){
      for(int e1 = 0; e1 < ne; e1++){
        evol_b1_proj[rho][flv][e1] = b0_proj[flv];

        AntineutrinoCPFix(rho);
        evol_b1_proj[rho][flv][e1].RotateToB1(params);
        AntineutrinoCPFix(rho);
      }
      b1_proj[rho][flv] = b0_proj[flv];

      AntineutrinoCPFix(rho);
      b1_proj[rho][flv].RotateToB1(params);
      AntineutrinoCPFix(rho);
    }
  }
}

SU_vector nuSQUIDS::GetState(int ie, int rho){
  return state[ie].rho[rho];
}

SU_vector nuSQUIDS::GetFlavorProj(int flv,int rho){
  return b1_proj[rho][flv];
}

SU_vector nuSQUIDS::GetMassProj(int flv,int rho){
  return b0_proj[flv];
}

SU_vector nuSQUIDS::GetHamiltonian(std::shared_ptr<Track> track, double E, int rho){
  index_rho = rho;
  EvolveProjectors(track->GetX());
  return H0(E)+HI(track->GetX(),E);
}

void nuSQUIDS::WriteStateHDF5(std::string str,std::string grp){

  hid_t error_stack;
  //H5Eset_auto(error_stack, NULL, NULL);
  H5Eset_auto (H5E_DEFAULT,NULL, NULL);

  hid_t file_id,group_id,root_id;
  hid_t dset_id;
  herr_t  status;
  // create HDF5 file
  //std::cout << "writing to hdf5 file" << std::endl;
  // H5F_ACC_TRUNC : overwrittes file
  // H5F_ACC_EXCL  : files if file existsi 
  file_id = H5Fopen(str.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0 ) // file already exists
    file_id = H5Fcreate(str.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  root_id = H5Gopen(file_id, "/",H5P_DEFAULT);
  if ( grp != "/" )
    group_id = H5Gcreate(root_id, grp.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  else
    group_id = root_id;

  // write the energy range
  hsize_t Edims[1]={E_range.size()};
  dset_id = H5LTmake_dataset(group_id,"energies",1,Edims,H5T_NATIVE_DOUBLE,E_range.data());
  H5LTset_attribute_string(group_id, "energies", "elogscale", (elogscale) ? "True":"False");

  // write mixing parameters
  hsize_t dim[1]{1};
  H5LTmake_dataset(group_id,"basic",1,dim,H5T_NATIVE_DOUBLE,0);
  H5LTmake_dataset(group_id,"mixingangles",1,dim,H5T_NATIVE_DOUBLE,0);
  H5LTmake_dataset(group_id,"CPphases",1,dim,H5T_NATIVE_DOUBLE,0);
  H5LTmake_dataset(group_id,"massdifferences",1,dim,H5T_NATIVE_DOUBLE,0);

  H5LTset_attribute_int(group_id, "basic","numneu",&numneu, 1);
  H5LTset_attribute_string(group_id, "basic","NT",NT.c_str());
  H5LTset_attribute_string(group_id, "basic", "interactions", (iinteraction) ? "True":"False");

  for ( int i = 0; i < 15; i++){
    std::string label = param_label_map[i];
    double value = params.GetMixingAngle(param_label_index[i][0],param_label_index[i][1]);
    //gsl_matrix_get(params.th, param_label_index[i][0],param_label_index[i][1]);
    H5LTset_attribute_double(group_id, "mixingangles",label.c_str(),&value, 1);
  }
  for ( int i = 15; i < 18; i++){
    std::string label = param_label_map[i];
    double value = params.GetPhase(param_label_index[i][0],param_label_index[i][1]);
    //gsl_matrix_get(params.dcp, param_label_index[i][0],param_label_index[i][1]);
    H5LTset_attribute_double(group_id, "CPphases",label.c_str(),&value, 1);
  }
  for ( int i = 18; i < 23; i++){
    std::string label = param_label_map[i];
    double value = params.GetSquaredEnergyDifference(param_label_index[i][0]);
    //gsl_matrix_get(params.dmsq, param_label_index[i][0],param_label_index[i][1]);
    H5LTset_attribute_double(group_id, "massdifferences",label.c_str(),&value, 1);
  }
  //writing state
  const int numneusq = numneu*numneu;
  hsize_t statedim[2] {E_range.size(),(hsize_t)numneu*numneu};
  std::vector<double> neustate(numneusq*ne), aneustate(numneusq*ne);

  for(int ie = 0; ie < ne; ie++){
    for(int i = 0; i < numneu*numneu; i ++){
        if (NT == "both"){
          neustate[ie*numneusq + i] = state[ie].rho[0][i];
          aneustate[ie*numneusq + i] = state[ie].rho[1][i];
        }
        else if (NT == "neutrino"){
          neustate[ie*numneusq + i] = state[ie].rho[0][i];
          aneustate[ie*numneusq + i] = 0.0;
        }
        else if (NT == "antineutrino"){
          neustate[ie*numneusq + i] = 0.0;
          aneustate[ie*numneusq + i] = state[ie].rho[0][i];
        }
      }
    }

  dset_id = H5LTmake_dataset(group_id,"neustate",2,statedim,H5T_NATIVE_DOUBLE,(void*)neustate.data());
  dset_id = H5LTmake_dataset(group_id,"aneustate",2,statedim,H5T_NATIVE_DOUBLE,(void*)aneustate.data());

  // writing state flavor and mass composition
  hsize_t pdim[2] {E_range.size(), (hsize_t) numneu};
  if ( NT == "both" )
    pdim[1] *= pdim[1];
  std::vector<double> flavor,mass;

  for(int ie = 0; ie < ne; ie++){
    // neutrino
    for(int i = 0; i < numneu; i++){
      if (NT == "both" or NT == "neutrino"){
          flavor.push_back(EvalFlavorAtNode(i,ie,0));
          mass.push_back(EvalMassAtNode(i,ie,0));
        }
        else if (NT == "antineutrino"){
          neustate.push_back(0.0);
          aneustate.push_back(0.0);
        }
    }
      // antineutrino
    for(int i = 0; i < numneu; i++){
      if (NT == "both" or NT == "antineutrino"){
          flavor.push_back(EvalFlavorAtNode(i,ie,0));
          mass.push_back(EvalMassAtNode(i,ie,0));
        }
        else if (NT == "neutrino"){
          neustate.push_back(0.0);
          aneustate.push_back(0.0);
        }
    }
  }
  dset_id = H5LTmake_dataset(group_id,"flavorcomp",2,pdim,H5T_NATIVE_DOUBLE,(void*)flavor.data());
  dset_id = H5LTmake_dataset(group_id,"masscomp",2,pdim,H5T_NATIVE_DOUBLE,(void*)mass.data());

  // writing body and track information
  hsize_t trackparamdim[1] {track->GetTrackParams().size()};
  if ( trackparamdim[0] == 0 ) {
    H5LTmake_dataset(group_id,"track",1,dim,H5T_NATIVE_DOUBLE,0);
  } else {
    H5LTmake_dataset(group_id,"track",1,trackparamdim,H5T_NATIVE_DOUBLE,track->GetTrackParams().data());
  }

  double xi = track->GetInitialX();
  H5LTset_attribute_double(group_id, "track","XINI",&xi, 1);
  double xf = track->GetFinalX();
  H5LTset_attribute_double(group_id, "track","XEND",&xf, 1);
  double xx = track->GetX();
  H5LTset_attribute_double(group_id, "track","X",&xx, 1);

  hsize_t bodyparamdim[1] {body->GetBodyParams().size()};
  if ( bodyparamdim[0] == 0 ){
    H5LTmake_dataset(group_id,"body",1,dim,H5T_NATIVE_DOUBLE,0);
  } else {
    H5LTmake_dataset(group_id,"body",1,bodyparamdim,H5T_NATIVE_DOUBLE,body->GetBodyParams().data());
  }
  H5LTset_attribute_string(group_id, "body", "NAME", body->name.c_str());
  int bid = body->id;
  H5LTset_attribute_int(group_id, "body", "ID", &bid,1);

  // writing cross section information 
  //
  //
  // close HDF5 file
  H5Gclose ( root_id );
  if ( root_id != group_id )
    H5Gclose ( group_id );
  H5Fclose (file_id);

}

void nuSQUIDS::ReadStateHDF5(std::string str,std::string grp){
  hid_t file_id,group_id,root_id,status;
  // open HDF5 file
  //std::cout << "reading from hdf5 file" << std::endl;
  file_id = H5Fopen(str.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  group_id = H5Gopen(root_id, grp.c_str(), H5P_DEFAULT);
  if ( group_id < 0 )
      throw std::runtime_error("nuSQUIDS::Error::Group '" + grp + "' does not exist in HDF5.");

  // read number of neutrinos
  H5LTget_attribute_int(group_id, "basic", "numneu", &numneu);
  // neutrino/antineutrino/both
  char auxchar[20];
  H5LTget_attribute_string(group_id, "basic", "NT", auxchar);
  NT = auxchar;
  // interactions
  H5LTget_attribute_string(group_id,"basic","interactions", auxchar);
  std::string aux = auxchar;
  if ( aux == "True")
    iinteraction = true;
  else
    iinteraction = false;

  // read and set mixing parameters
  for ( int i = 0 ; i < 15; i++){
    double value;
    H5LTget_attribute_double(group_id,"mixingangles", param_label_map[i].c_str(), &value);
    //Set(param_label_map[i], value);
    Set(static_cast<MixingParameter>(i), value);
  }
  for ( int i = 15 ; i < 18; i++){
    double value;
    H5LTget_attribute_double(group_id,"CPphases", param_label_map[i].c_str(), &value);
    //Set(param_label_map[i], value);
    Set(static_cast<MixingParameter>(i), value);
  }
  for ( int i = 18 ; i < 23; i++){
    double value;
    H5LTget_attribute_double(group_id,"massdifferences", param_label_map[i].c_str(), &value);
    //Set(param_label_map[i], value);
    Set(static_cast<MixingParameter>(i), value);
  }

  // reading energy
  hsize_t dims[2];
  H5LTget_dataset_info(group_id, "energies", dims, NULL, NULL);

  double data[dims[0]];
  ne = dims[0];
  H5LTread_dataset_double(group_id, "energies", data);
  //for (int i = 0; i < dims[0]; i ++ )
  //  std::cout << data[i] << std::endl;

  H5LTget_attribute_string(group_id,"energies","elogscale", auxchar);
  aux = auxchar;
  if ( aux == "True")
    elogscale = true;
  else
    elogscale = false;

  // initializing nuSQUIDS
  if (ne == 1){
    if(not inusquids)
      init();
    Set_E(data[0]);
  }
  else {
    init(data[0]/units.GeV,data[ne-1]/units.GeV,ne);
  }

  // reading state
  H5LTget_dataset_info(group_id,"neustate", dims,NULL,NULL);
  double neudata[dims[0]*dims[1]];
  H5LTread_dataset_double(group_id,"neustate", neudata);

  H5LTget_dataset_info(group_id,"aneustate", dims,NULL,NULL);
  double aneudata[dims[0]*dims[1]];
  H5LTread_dataset_double(group_id,"aneustate", aneudata);

  for(int ie = 0; ie < dims[0]; ie++){
    for (int j = 0; j < dims[1]; j ++){
      if (NT == "neutrino")
        state[ie].rho[0][j] = neudata[ie*dims[1]+j];
      else if ( NT == "antineutrino")
        state[ie].rho[0][j] = aneudata[ie*dims[1]+j];
      else if ( NT == "both" ){
        state[ie].rho[0][j] = neudata[ie*dims[1]+j];
        state[ie].rho[1][j] = aneudata[ie*dims[1]+j];
      }
    }
  }

  // reading body and track
  int id;
  hsize_t dimbody[2];
  H5LTget_attribute_int(group_id,"body","ID",&id);

  H5LTget_dataset_info(group_id,"body", dimbody,NULL,NULL);
  double body_params[dimbody[0]];
  H5LTread_dataset_double(group_id,"body", body_params);

  hsize_t dimtrack[2];
  H5LTget_dataset_info(group_id,"track", dimtrack ,NULL,NULL);
  double track_params[dimtrack[0]];
  H5LTread_dataset_double(group_id,"track", track_params);
  double x_current;
  H5LTget_attribute_double(group_id,"track","X",&x_current);


  // setting body and track
  SetBodyTrack(id,dimbody[0],body_params,dimtrack[0],track_params);

  // set projector to current position
  track->SetX(x_current);
  EvolveProjectors(track->GetX());
  t = track->GetX();
  t_ini = track->GetInitialX();

  // close HDF5 file
  H5Gclose ( group_id );
  H5Gclose ( root_id );
  H5Fclose (file_id);
}


void nuSQUIDS::SetBodyTrack(int body_id, int body_params_len, double body_params[], int track_params_len, double track_params[]){
    switch(body_id){
      case 1:
        {
          body = std::make_shared<Vacuum>();
          track = std::make_shared<Vacuum::Track>(track_params[0],track_params[1]);
          break;
        }
      case 2:
        {
          body = std::make_shared<ConstantDensity>(body_params[0],body_params[1]);
          track = std::make_shared<ConstantDensity::Track>(track_params[0],track_params[1]);
          break;
        }
      case 3:
        {
          const int xn = body_params_len/3;
          std::vector<double> xx(xn),rho(xn),ye(xn);
          for(int i = 0; i < xn; i++){
            xx[i] = body_params[i];
            rho[i] = body_params[xn+i];
            ye[i] = body_params[2*xn+i];
          }
          body = std::make_shared<VariableDensity>(xx,rho,ye);
          track = std::make_shared<VariableDensity::Track>(track_params[0],track_params[1]);
          break;
        }
      case 4:
        {
          body = std::make_shared<Earth>();
          track = std::make_shared<Earth::Track>(track_params[0],track_params[1],track_params[2]);
          break;
        }
      case 5:
        {
          body = std::make_shared<Sun>();
          track = std::make_shared<Sun::Track>(track_params[0],track_params[1]);
          break;
        }
      case 6:
        {
          body = std::make_shared<SunASnu>();
          track = std::make_shared<SunASnu::Track>(track_params[0],track_params[1]);
          break;
        }
      case 7:
        {
          body = std::make_shared<EarthAtm>();
          track = std::make_shared<EarthAtm::Track>(track_params[0]);
          break;
        }
      default:
        {
          std::cerr << "nuSQUIDS::SetBodyTrack : unknown body/track" << std::endl;
          exit(1);
          break;
        }
    }
}

int nuSQUIDS::GetNumNeu(){
  return numneu;
}

void nuSQUIDS::ProgressBar(){
  double progress = track->GetX()/track->GetFinalX();
  int barWidth = 70;
  int pos = barWidth * progress;
  std::cout << "[";
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();

}

void nuSQUIDS::Set_TauRegeneration(bool opt){
    if ( NT != "both" and opt )
      throw std::runtime_error("nuSQUIDS::Error::Cannot set TauRegeneration to True when NT != 'both'.");
    tauregeneration = opt;
}

void nuSQUIDS::Set_ProgressBar(bool opt){
    progressbar = opt;
}

std::shared_ptr<Track> nuSQUIDS::GetTrack(void){
  return track;
}

std::shared_ptr<Body> nuSQUIDS::GetBody(void){
  return body;
}

void nuSQUIDS::Set(MixingParameter p, double val){
  switch (p) {
    case TH12:
      params.SetMixingAngle(0,1,val);
      break;
    case TH13:
      params.SetMixingAngle(0,2,val);
      break;
    case TH23:
      params.SetMixingAngle(1,2,val);
      break;
    case TH14:
      params.SetMixingAngle(0,3,val);
      break;
    case TH24:
      params.SetMixingAngle(1,3,val);
      break;
    case TH34:
      params.SetMixingAngle(2,3,val);
      break;
    case TH15:
      params.SetMixingAngle(0,4,val);
      break;
    case TH25:
      params.SetMixingAngle(1,4,val);
      break;
    case TH35:
      params.SetMixingAngle(2,4,val);
      break;
    case TH45:
      params.SetMixingAngle(3,4,val);
      break;
    case TH16:
      params.SetMixingAngle(0,5,val);
      break;
    case TH26:
      params.SetMixingAngle(1,5,val);
      break;
    case TH36:
      params.SetMixingAngle(2,5,val);
      break;
    case TH46:
      params.SetMixingAngle(3,5,val);
      break;
    case TH56:
      params.SetMixingAngle(4,5,val);
      break;
    case DM21SQ:
      params.SetSquaredEnergyDifference(1,val);
      break;
    case DM31SQ:
      params.SetSquaredEnergyDifference(2,val);
      break;
    case DM41SQ:
      params.SetSquaredEnergyDifference(3,val);
      break;
    case DM51SQ:
      params.SetSquaredEnergyDifference(4,val);
      break;
    case DM61SQ:
      params.SetSquaredEnergyDifference(5,val);
      break;
    case DELTA1:
      params.SetPhase(0,2,val);
      break;
    case DELTA2:
      params.SetPhase(0,4,val);
      break;
    case DELTA3:
      params.SetPhase(0,5,val);
      break;
  }
}

void nuSQUIDS::Set_MixingParametersToDefault(void){
  // set parameters as in arXiv:1409.5439 NO
  // but with delta_CP = 0.0
  Set(TH12,0.583996);
  Set(TH13,0.148190);
  Set(TH23,0.737324);

  Set(DM21SQ,7.5e-05);
  Set(DM31SQ,0.00257);

  Set(DELTA1,0.0);
}

/*
nuSQUIDS::~nuSQUIDS(void){
  if (inusquids){
    free();
  }
}
*/

void nuSQUIDS::Set_Basis(BASIS b){
  basis = b;
}


//==================================================================
//==================================================================
//==================================================================

nuSQUIDSAtm::nuSQUIDSAtm(double costh_min,double costh_max,int costh_div,
                         double energy_min,double energy_max,int energy_div,
                         int numneu,std::string NT,
                         bool elogscale,bool iinteraction){

  nusq_array = std::vector<nuSQUIDS>(costh_div);
  costh_array = linspace(costh_min,costh_max,costh_div-1);
  if(elogscale){
    enu_array = logspace(energy_min,energy_max,energy_div-1);
  }
  else{
    enu_array = linspace(energy_min,energy_max,energy_div-1);
  }
  std::transform(enu_array.begin(),enu_array.end(),std::back_inserter(log_enu_array),[](int enu){ return log(enu);});

  earth_atm = std::make_shared<EarthAtm>();
  for(double costh : costh_array)
    track_array.push_back(std::make_shared<EarthAtm::Track>(acos(costh)));

  int i = 0;
  for(nuSQUIDS& nsq : nusq_array){
    nsq.Init(energy_min,energy_max,energy_div,numneu,NT,elogscale,interaction);
    nsq.Set_Body(earth_atm);
    nsq.Set_Track(track_array[i]);
    i++;
  }

  inusquidsatm = true;
}

void nuSQUIDSAtm::EvolveState(void){
  if(not iinistate)
    throw std::runtime_error("nuSQUIDSAtm::Error::State not initialized.");
  if(not inusquidsatm)
    throw std::runtime_error("nuSQUIDSAtm::Error::nuSQUIDSAtm not initialized.");
  int i = 0;
  for(nuSQUIDS& nsq : nusq_array){
    if(progressbar){
      std::cout << "Calculating cos(th) = " + std::to_string(costh_array[i]) << std::endl;
      i++;
    }
    nsq.EvolveState();
    if(progressbar){
      std::cout << std::endl;
    }
  }
}

void nuSQUIDSAtm::Set_initial_state(array3D ini_flux, std::string basis){
  if(ini_flux.size() != costh_array.size())
    throw std::runtime_error("nuSQUIDSAtm::Error::First dimension of input array is incorrect.");
  int i = 0;
  for(nuSQUIDS& nsq : nusq_array){
    nsq.Set_initial_state(ini_flux[i],basis);
    i++;
  }
  iinistate = true;
}

void nuSQUIDSAtm::Set_initial_state(array4D ini_flux, std::string basis){
  if(ini_flux.size() != costh_array.size())
    throw std::runtime_error("nuSQUIDSAtm::Error::First dimension of input array is incorrect.");
  int i = 0;
  for(nuSQUIDS& nsq : nusq_array){
    nsq.Set_initial_state(ini_flux[i],basis);
    i++;
  }
  iinistate = true;
}

void nuSQUIDSAtm::WriteStateHDF5(std::string filename){
  if(not iinistate)
    throw std::runtime_error("nuSQUIDSAtm::Error::State not initialized.");
  if(not inusquidsatm)
    throw std::runtime_error("nuSQUIDSAtm::Error::nuSQUIDSAtm not initialized.");

  hid_t file_id,group_id,root_id;
  hid_t dset_id;
  // create HDF5 file
  file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  root_id = H5Gopen(file_id, "/",H5P_DEFAULT);

  // write the zenith range
  hsize_t costhdims[1]={costh_array.size()};
  dset_id = H5LTmake_dataset(root_id,"zenith_angles",1,costhdims,H5T_NATIVE_DOUBLE,costh_array.data());
  hsize_t energydims[1]={enu_array.size()};
  dset_id = H5LTmake_dataset(root_id,"energy_range",1,energydims,H5T_NATIVE_DOUBLE,enu_array.data());

  H5Gclose (root_id);
  H5Fclose (file_id);

  int i = 0;
  for(nuSQUIDS& nsq : nusq_array){
    nsq.WriteStateHDF5(filename,"costh_"+std::to_string(costh_array[i]));
    i++;
  }
}


void nuSQUIDSAtm::Set_MixingParametersToDefault(){
  for(nuSQUIDS& nsq : nusq_array){
    nsq.Set_MixingParametersToDefault();
  }
}

void nuSQUIDSAtm::Set(MixingParameter p,double v){
  for(nuSQUIDS& nsq : nusq_array){
    nsq.Set(p,v);
  }
}

void nuSQUIDSAtm::Set_TauRegeneration(bool v){
  for(nuSQUIDS& nsq : nusq_array){
    nsq.Set_TauRegeneration(v);
  }
}

void nuSQUIDSAtm::ReadStateHDF5(std::string filename){

  hid_t file_id,group_id,root_id;
  hid_t dset_id;
  // create HDF5 file
  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  root_id = H5Gopen(file_id, "/",H5P_DEFAULT);
  group_id = root_id;

  // read the zenith range dimension
  hsize_t costhdims[1];
  H5LTget_dataset_info(group_id, "zenith_angles", costhdims, NULL, NULL);

  double data[costhdims[0]];
  H5LTread_dataset_double(group_id, "zenith_angles", data);
  costh_array.clear();
  costh_array.resize(costhdims[0]);
  for (int i = 0; i < costhdims[0]; i ++)
    costh_array[i] = data[i];

  hsize_t energydims[1];
  H5LTget_dataset_info(group_id, "energy_range", energydims, NULL, NULL);

  double enu_data[energydims[0]];
  H5LTread_dataset_double(group_id, "energy_range", enu_data);
  enu_array.clear();log_enu_array.clear();
  enu_array.resize(energydims[0]);log_enu_array.resize(energydims[0]);
  for (int i = 0; i < energydims[0]; i ++){
    enu_array[i] = enu_data[i];
    log_enu_array[i] = log(enu_data[i]);
  }

  H5Gclose(root_id);
  H5Fclose(file_id);

  // resize apropiately the nuSQUIDSAtm container vector
  nusq_array.clear();
  nusq_array = std::vector<nuSQUIDS>(costhdims[0]);

  int i = 0;
  for(nuSQUIDS& nsq : nusq_array){
    nsq.ReadStateHDF5(filename,"costh_"+std::to_string(costh_array[i]));
    i++;
  }

  iinistate = true;
  inusquidsatm = true;
}

double nuSQUIDSAtm::LinInter(double x,double xM, double xP, double yM, double yP) const {
  return yM + (yP-yM)*(x-xM)/(xP-xM);
}

double nuSQUIDSAtm::EvalFlavor(int flv,double costh,double enu,int rho){
  // here the energy enters in GeV
  if(not iinistate)
    throw std::runtime_error("nuSQUIDSAtm::Error::State not initialized.");
  if(not inusquidsatm)
    throw std::runtime_error("nuSQUIDSAtm::Error::nuSQUIDSAtm not initialized.");

  if( costh < costh_array[0] or costh > costh_array.back())
    throw std::runtime_error("nuSQUIDSAtm::Error::EvalFlavor::cos(th) out of bounds.");
  if( enu < enu_array[0] or enu > enu_array.back() )
    throw std::runtime_error("nuSQUIDSAtm::Error::EvalFlavor::neutrino energy out of bounds.");

  std::shared_ptr<EarthAtm::Track> track = std::make_shared<EarthAtm::Track>(acos(costh));
  // get the evolution generator
  SU_vector H0_at_enu = nusq_array[0].H0(enu*units.GeV);
  // get the evolved projector for the right distance and energy
  SU_vector evol_proj = nusq_array[0].GetFlavorProj(flv,rho).SUEvolve(H0_at_enu,track->GetFinalX());

  int cth_M = -1;
  for(int i = 0; i < costh_array.size(); i++){
    if ( costh >= costh_array[i] and costh <= costh_array[i+1] ) {
      cth_M = i;
      break;
    }
  }

  int loge_M = -1;
  double logE = log(enu);
  for(int i = 0; i < log_enu_array.size(); i++){
    if ( logE >= log_enu_array[i] and logE <= log_enu_array[i+1] ) {
      loge_M = i;
      break;
    }
  }

  //std::cout << cth_M << " " << loge_M << std::endl;

  double phiMM,phiMP,phiPM,phiPP;
  phiMM = nusq_array[cth_M].GetState(loge_M,rho)*evol_proj;
  phiMP = nusq_array[cth_M].GetState(loge_M+1,rho)*evol_proj;
  phiPM = nusq_array[cth_M+1].GetState(loge_M,rho)*evol_proj;
  phiPP = nusq_array[cth_M+1].GetState(loge_M+1,rho)*evol_proj;

  return LinInter((double)costh,(double)costh_array[cth_M],(double)costh_array[cth_M+1],
        (double)LinInter((double)logE,(double)log_enu_array[loge_M],(double)log_enu_array[loge_M+1],(double)phiMM,(double)phiMP),
        (double)LinInter((double)logE,(double)log_enu_array[loge_M],(double)log_enu_array[loge_M+1],(double)phiPM,(double)phiPP));
}

void nuSQUIDSAtm::Set_rel_error(double er){
  for(nuSQUIDS& nsq : nusq_array){
    nsq.Set_rel_error(er);
  }
}

void nuSQUIDSAtm::Set_abs_error(double er){
  for(nuSQUIDS& nsq : nusq_array){
    nsq.Set_abs_error(er);
  }
}

size_t nuSQUIDSAtm::GetNumE(void){
  return enu_array.size();
}
size_t nuSQUIDSAtm::GetNumCos(void){
  return costh_array.size();
}

void nuSQUIDSAtm::Set_ProgressBar(bool v){
    progressbar = v;
    for(nuSQUIDS& nsq : nusq_array){
      nsq.Set_ProgressBar(v);
    }
}

array1D nuSQUIDSAtm::GetERange(void){
  return enu_array;
}

array1D nuSQUIDSAtm::GetCosthRange(void){
  return costh_array;
}

} // close namespace
