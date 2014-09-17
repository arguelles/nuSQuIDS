#include "nuSQUIDS.h"

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
  ini(ne,numneu,1,0);
  Set("H1",true);
  Set("h_max", std::numeric_limits<double>::max() );

  //===============================
  // init projectors             //
  //===============================

  b0_proj.resize(numneu);
  for(int flv = 0; flv < numneu; flv++){
    b0_proj[flv].InitSU_vector("Proj",flv,nsun);
  }

  b1_proj.resize(nrhos);
  for(int rho = 0; rho < nrhos; rho++){
    b1_proj[rho].resize(numneu);
    for(int flv = 0; flv < numneu; flv++){
      b1_proj[rho][flv].InitSU_vector("Proj",flv,nsun);
      b1_proj[rho][flv].RotateToB1(&params);
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
        evol_b0_proj[rho][flv][e1].InitSU_vector("Proj",flv,nsun);
        evol_b1_proj[rho][flv][e1].InitSU_vector("Proj",flv,nsun);
        evol_b1_proj[rho][flv][e1].RotateToB1(&params);
      }
    }
  }

  //===============================
  // init square mass difference //
  //===============================

  iniH0();

  //===============================
  // END                         //
  //===============================
  inusquids = true;
}

void nuSQUIDS::Set_E(double Enu){
  if( ne != 1 )
    throw std::runtime_error("nuSQUIDS::Error:Cannot use Set_E in single energy mode.");
  E_range = vector<double>{Enu};
  Set_xrange(E_range[0],E_range[ne-1],"lin");

  ienergy = true;
}

/*
void nuSQUIDS::Set_E(vector<double> Enu){

}
*/

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
    ini(ne,numneu,nrhos,nrhos);
  else
    ini(ne,numneu,nrhos,0);

  SetScalarsToZero();

  t = 0;

  Set("H1",true);
  Set("h_max", std::numeric_limits<double>::max() );

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
  // init projectors             //
  //===============================

  b0_proj.resize(numneu);
  for(int flv = 0; flv < numneu; flv++){
    b0_proj[flv].InitSU_vector("Proj",flv,numneu);
  }

  b1_proj.resize(nrhos);
  for(int rho = 0; rho < nrhos; rho++){
    b1_proj[rho].resize(numneu);
    for(int flv = 0; flv < numneu; flv++){
      b1_proj[rho][flv].InitSU_vector("Proj",flv,nsun);
      b1_proj[rho][flv].RotateToB1(&params);
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
        evol_b0_proj[rho][flv][e1].InitSU_vector("Proj",flv,nsun);
        evol_b1_proj[rho][flv][e1].InitSU_vector("Proj",flv,nsun);
        evol_b1_proj[rho][flv][e1].RotateToB1(&params);
      }
    }
  }

  //===============================
  // init square mass difference //
  //===============================

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
    Set("NonCohInt",true);
    Set("ScalInt",true);
    Set("OInt",true);
  }

  //===============================
  // END                         //
  //===============================

  inusquids = true;
}

void nuSQUIDS::PreDerive(double x){
  track->SetX(tunit*x);
  EvolveProjectors(tunit*x);
  if(iinteraction){
    UpdateInteractions();
  }
}

void nuSQUIDS::EvolveProjectors(double x){
  for(int rho = 0; rho < nrhos; rho++){
    for(int flv = 0; flv < numneu; flv++){
      for(int ei = 0; ei < ne; ei++){
        SU_vector h0 = H0(E_range[ei]);
        evol_b0_proj[rho][flv][ei] = b0_proj[flv].SUEvolve(h0,tunit*(x-t_ini));
        evol_b1_proj[rho][flv][ei] = b1_proj[rho][flv].SUEvolve(h0,tunit*(x-t_ini));
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
    SU_vector potential = (CC+NC)*evol_b1_proj[index_rho][0][ei] + (NC)*(evol_b1_proj[index_rho][1][ei] + evol_b1_proj[index_rho][2][ei]);

    if (index_rho == 0){
        // neutrino potential
        return potential;
    } else if (index_rho == 1){
        // antineutrino potential
        return (-1.0)*potential;
    } else{
        std::cerr << "nuSQUIDS::HI : unknown particle or antiparticle" << std::endl;
        exit(1);
    }
}

SU_vector nuSQUIDS::GammaRho(int ei){
    if (not iinteraction){
    SU_vector V(nsun);
    return V;
    }

    //std::cout << invlen_INT[index_rho][0][ei] << " " << invlen_INT[index_rho][0][ei] << " " << invlen_INT[index_rho][0][ei] << std::endl;
    //std::cout << "GammaRho" << std::endl;
    //std::cout << evol_b1_proj[index_rho][0][ei]*(0.5*invlen_INT[index_rho][0][ei]) +
    //       evol_b1_proj[index_rho][1][ei]*(0.5*invlen_INT[index_rho][1][ei]) +
    //       evol_b1_proj[index_rho][2][ei]*(0.5*invlen_INT[index_rho][2][ei]) << std::endl;

    return evol_b1_proj[index_rho][0][ei]*(0.5*invlen_INT[index_rho][0][ei]) +
           evol_b1_proj[index_rho][1][ei]*(0.5*invlen_INT[index_rho][1][ei]) +
           evol_b1_proj[index_rho][2][ei]*(0.5*invlen_INT[index_rho][2][ei]);
}

SU_vector nuSQUIDS::InteractionsRho(int e1){
  if (not iinteraction){
    SU_vector V(nsun);
    return V;
  }

  // this implements the NC interactinos
  // the tau regeneration terms are implemented at the end
  SU_vector nc_term(nsun);
  for(int e2 = e1 + 1; e2 < ne; e2++){
    // here we assume the cross section to be the same for all flavors
    //std::cout << dNdE_NC[index_rho][0][e2][e1] << " " << invlen_NC[index_rho][0][e2] << std::endl;
    nc_term += SU.ACommutator(evol_b1_proj[index_rho][0][e1] + evol_b1_proj[index_rho][1][e1] + evol_b1_proj[index_rho][2][e1],
                              state[e2].rho[index_rho])*
               (0.5*dNdE_NC[index_rho][0][e2][e1]*invlen_NC[index_rho][0][e2]);
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
    EvolveSUN(track->GetInitialX(),track->GetFinalX());
  }
  else {
    int tau_steps = (int) ((track->GetFinalX() - track->GetInitialX())/tau_reg_scale);
    std::cout << tau_steps << " " << tau_reg_scale/units.km<< std::endl;
    for (int i = 0; i < tau_steps; i++){
      double x_inter = track->GetInitialX() + (double)(i)*tau_reg_scale;
      std::cout << x_inter/units.km << std::endl;
      EvolveSUN(x_inter,x_inter + tau_reg_scale);
      ConvertTauIntoNuTau();
    }
    EvolveSUN(tau_reg_scale*tau_steps,track->GetFinalX());
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
  t_end = 0.0;
  t_ini = 0.0;
  t = 0.0;
}

void nuSQUIDS::Set_initial_state(vector<double> v, string basis){
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

void nuSQUIDS::Set_initial_state(array2D v, string basis){
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

void nuSQUIDS::Set_initial_state(array3D v, string basis){
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

void nuSQUIDS::WriteState(string filename){
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

void nuSQUIDS::ReadState(string filename){

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
  return GetExpectationValueD(b0_proj[flv], rho, EE);
}

double nuSQUIDS::EvalFlavor(int flv,double EE, int rho){
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != "both" )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  return GetExpectationValueD(b1_proj[rho][flv], rho, EE);
}

double nuSQUIDS::EvalMassAtNode(int flv, int ei, int rho){
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != "both" )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  return GetExpectationValue(b0_proj[flv], rho, ei);
}

double nuSQUIDS::EvalFlavorAtNode(int flv, int ei, int rho){
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != "both" )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
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

  return GetExpectationValue(b1_proj[0][flv], 0, 0);
}

void nuSQUIDS::iniH0(){
  DM2.InitSU_vector(nsun);
  for(int i = 1; i < nsun; i++){
      DM2 += (b0_proj[i])*gsl_matrix_get(params.dmsq,i,0);
  }
}

void nuSQUIDS::iniProyectors(){

  for(int rho = 0; rho < nrhos; rho++){
    for(int flv = 0; flv < numneu; flv++){
      b1_proj[rho][flv].InitSU_vector("Proj",flv,nsun);
      b1_proj[rho][flv].RotateToB1(&params);
    }
  }

  for(int rho = 0; rho < nrhos; rho++){
    for(int flv = 0; flv < numneu; flv++){
      for(int e1 = 0; e1 < ne; e1++){
        evol_b1_proj[rho][flv][e1].InitSU_vector("Proj",flv,nsun);
        evol_b1_proj[rho][flv][e1].RotateToB1(&params);
      }
    }
  }

}

void nuSQUIDS::SetIniFlavorProyectors(){
  for(int rho = 0; rho < nrhos; rho++){
    for(int flv = 0; flv < numneu; flv++){
      for(int e1 = 0; e1 < ne; e1++){
        evol_b1_proj[rho][flv][e1] = b0_proj[flv];
        evol_b1_proj[rho][flv][e1].RotateToB1(&params);
      }
      b1_proj[rho][flv] = b0_proj[flv];
      b1_proj[rho][flv].RotateToB1(&params);
    }
  }
}

SU_vector nuSQUIDS::GetHamiltonian(std::shared_ptr<Track> track, double E, int rho){
  index_rho = rho;
  EvolveProjectors(track->GetX());
  return H0(E)+HI(track->GetX(),E);
}

void nuSQUIDS::WriteStateHDF5(string str){

  hid_t file_id;
  hid_t   dset_id;
  herr_t  status;
  // create HDF5 file
  file_id = H5Fcreate (str.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // write the energy range
  hsize_t Edims[1]={E_range.size()};
  dset_id = H5LTmake_dataset(file_id,"/energies",1,Edims,H5T_NATIVE_DOUBLE,E_range.data());
  H5LTset_attribute_string(file_id, "/energies", "elogscale", (elogscale) ? "True":"False");


  // write mixing parameters
  hsize_t dim[1]{1};
  H5LTmake_dataset(file_id,"/basic",1,dim,H5T_NATIVE_DOUBLE,0);
  H5LTmake_dataset(file_id,"/mixingangles",1,dim,H5T_NATIVE_DOUBLE,0);
  H5LTmake_dataset(file_id,"/CPphases",1,dim,H5T_NATIVE_DOUBLE,0);
  H5LTmake_dataset(file_id,"/massdifferences",1,dim,H5T_NATIVE_DOUBLE,0);

  H5LTset_attribute_int(file_id, "/basic","numneu",&numneu, 1);
  H5LTset_attribute_string(file_id, "/basic","NT",NT.c_str());
  H5LTset_attribute_string(file_id, "/basic", "interactions", (iinteraction) ? "True":"False");

  for ( int i = 0; i < 15; i++){
    string label = param_label_map[i];
    double value = gsl_matrix_get(params.th, param_label_index[i][0],param_label_index[i][1]);
    H5LTset_attribute_double(file_id, "/mixingangles",label.c_str(),&value, 1);
  }
  for ( int i = 15; i < 18; i++){
    string label = param_label_map[i];
    double value = gsl_matrix_get(params.dcp, param_label_index[i][0],param_label_index[i][1]);
    H5LTset_attribute_double(file_id, "/CPphases",label.c_str(),&value, 1);
  }
  for ( int i = 18; i < 23; i++){
    string label = param_label_map[i];
    double value = gsl_matrix_get(params.dmsq, param_label_index[i][0],param_label_index[i][1]);
    H5LTset_attribute_double(file_id, "/massdifferences",label.c_str(),&value, 1);
  }

  //writing state
  const int numneusq = numneu*numneu;
  hsize_t statedim[2] {E_range.size(),(hsize_t)numneu*numneu};
  vector<double> neustate, aneustate;

  for(int ie = 0; ie < ne; ie++){
    for(int i = 0; i < numneu*numneu; i ++){
        if (NT == "both"){
          neustate.push_back(state[ie].rho[0][i]);
          aneustate.push_back(state[ie].rho[1][i]);
        }
        else if (NT == "neutrino"){
          neustate.push_back(state[ie].rho[0][i]);
          aneustate.push_back(0.0);
        }
        else if (NT == "antineutrino"){
          neustate.push_back(0.0);
          aneustate.push_back(state[ie].rho[0][i]);
        }
      }
    }

  dset_id = H5LTmake_dataset(file_id,"/neustate",2,statedim,H5T_NATIVE_DOUBLE,neustate.data());
  dset_id = H5LTmake_dataset(file_id,"/aneustate",2,statedim,H5T_NATIVE_DOUBLE,aneustate.data());

  // writing state flavor and mass composition
  hsize_t pdim[2] {E_range.size(), (hsize_t) numneu};
  vector<double> flavor,mass;

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

  dset_id = H5LTmake_dataset(file_id,"/flavorcomp",2,pdim,H5T_NATIVE_DOUBLE,flavor.data());
  dset_id = H5LTmake_dataset(file_id,"/masscomp",2,pdim,H5T_NATIVE_DOUBLE,mass.data());

  // writing body and track information
  hsize_t trackparamdim[1] {track->GetTrackParams().size()};
  H5LTmake_dataset(file_id,"/track",1,trackparamdim,H5T_NATIVE_DOUBLE,track->GetTrackParams().data());

  double xx;
  xx = track->GetInitialX();
  H5LTset_attribute_double(file_id, "/track","XINI",&xx, 1);
  xx = track->GetFinalX();
  H5LTset_attribute_double(file_id, "/track","XEND",&xx, 1);
  xx = track->GetX();
  H5LTset_attribute_double(file_id, "/track","X",&xx, 1);

  hsize_t bodyparamdim[1] {body->GetBodyParams().size()};
  H5LTmake_dataset(file_id,"/body",1,bodyparamdim,H5T_NATIVE_DOUBLE,body->GetBodyParams().data());
  H5LTset_attribute_string(file_id, "/body", "NAME", body->name.c_str());
  int bid = body->id;
  H5LTset_attribute_int(file_id, "/body", "ID", &bid,1);

  // close HDF5 file
  status = H5Fclose (file_id);

}

void nuSQUIDS::ReadStateHDF5(string str){
  hid_t file_id, status;
  // open HDF5 file
  file_id = H5Fopen(str.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // read number of neutrinos
  H5LTget_attribute_int(file_id, "/basic", "numneu", &numneu);
  // neutrino/antineutrino/both
  char auxchar[20];
  H5LTget_attribute_string(file_id, "/basic", "NT", auxchar);
  NT = (string) auxchar;
  // interactions
  H5LTget_attribute_string(file_id,"/basic","interactions", auxchar);
  string aux = (string) auxchar;
  if ( aux == "True")
    iinteraction = true;
  else
    iinteraction = false;

  // read and set mixing parameters
  for ( int i = 0 ; i < 15; i++){
    double value;
    H5LTget_attribute_double(file_id,"/mixingangles", param_label_map[i].c_str(), &value);
    Set(param_label_map[i], value);
  }
  for ( int i = 15 ; i < 18; i++){
    double value;
    H5LTget_attribute_double(file_id,"/CPphases", param_label_map[i].c_str(), &value);
    Set(param_label_map[i], value);
  }
  for ( int i = 18 ; i < 23; i++){
    double value;
    H5LTget_attribute_double(file_id,"/massdifferences", param_label_map[i].c_str(), &value);
    Set(param_label_map[i], value);
  }

  // reading energy
  hsize_t dims[2];
  H5LTget_dataset_info(file_id, "/energies", dims, NULL, NULL);

  double data[dims[0]];
  ne = dims[0];
  H5LTread_dataset_double(file_id, "/energies", data);
  //for (int i = 0; i < dims[0]; i ++ )
  //  std::cout << data[i] << std::endl;

  H5LTget_attribute_string(file_id,"/energies","elogscale", auxchar);
  aux = (string) auxchar;
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
  H5LTget_dataset_info(file_id,"/neustate", dims,NULL,NULL);
  double neudata[dims[0]*dims[1]];
  H5LTread_dataset_double(file_id,"/neustate", neudata);

  H5LTget_dataset_info(file_id,"/aneustate", dims,NULL,NULL);
  double aneudata[dims[0]*dims[1]];
  H5LTread_dataset_double(file_id,"/aneustate", aneudata);

  for(int ie = 0; ie < dims[0]; ie++){
    for (int j = 0; j < dims[1]; j ++){
      if (NT == "neutrino")
        state[ie].rho[0][j] = neudata[ie*dims[1]+j];
      else if ( NT == "antineutrino")
        state[ie].rho[1][j] = aneudata[ie*dims[1]+j];
      else if ( NT == "both" ){
        state[ie].rho[0][j] = neudata[ie*dims[1]+j];
        state[ie].rho[1][j] = aneudata[ie*dims[1]+j];
      }
    }
  }

  // reading body and track
  int id;
  hsize_t dimbody[2];
  H5LTget_attribute_int(file_id,"/body","ID",&id);

  H5LTget_dataset_info(file_id,"/body", dimbody,NULL,NULL);
  double body_params[dimbody[0]];
  H5LTread_dataset_double(file_id,"/body", body_params);

  hsize_t dimtrack[2];
  H5LTget_dataset_info(file_id,"/track", dimtrack ,NULL,NULL);
  double track_params[dimtrack[0]];
  H5LTread_dataset_double(file_id,"/track", track_params);
  double x_current;
  H5LTget_attribute_double(file_id,"/track","X",&x_current);


  // setting body and track
  SetBodyTrack(id,dimbody[0],body_params,dimtrack[0],track_params);

  // set projector to current position
  track->SetX(x_current);
  EvolveProjectors(track->GetX());
  t = track->GetX();
  t_ini = track->GetInitialX();

  // close HDF5 file
  status = H5Fclose (file_id);
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
          vector<double> xx(xn),rho(xn),ye(xn);
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

void nuSQUIDS::Set_nuSQUIDS(string str, bool opt){
  if ( str == "tauregeneration" | str == "TauRegeneration" | str == "taureg" ) {
    if ( NT != "both" and opt )
      throw std::runtime_error("nuSQUIDS::Error::Cannot set TauRegeneration to True when NT != 'both'.");
    tauregeneration = opt;
  } else  {
    throw std::runtime_error("nuSQUIDS::Error::Unknown option.");
  }
}

std::shared_ptr<Track> nuSQUIDS::GetTrack(void){
  return track;
}

std::shared_ptr<Body> nuSQUIDS::GetBody(void){
  return body;
}

/*
nuSQUIDS::~nuSQUIDS(void){
  if (inusquids){
    free();
  }
}
*/
