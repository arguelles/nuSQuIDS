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

#define H5Gopen_vers 2
#define H5Gcreate_vers 2
#define H5Eset_auto_vers 2

#include "nuSQuIDS.h"

namespace nusquids{

void nuSQUIDS::init(double xini){
  // single energy implementation
  ne = 1;

  if (NT == neutrino || NT == antineutrino)
    nrhos = 1;
  else {
    throw std::runtime_error("nuSQUIDS::Error::NT = {neutrino,antineutrino} not : " + std::to_string(NT));
  }

  if ( numneu > SQUIDS_MAX_HILBERT_DIM )
    throw std::runtime_error("nuSQUIDS::Error::Maximum number of neutrinos exceded");
  nsun = numneu;

  //initialize SQUIDS
  ini(ne,numneu,1,0,xini);
  Set_CoherentRhoTerms(true);
  Set_h_max(std::numeric_limits<double>::max() );

  //===============================
  // set parameters to default   //
  //===============================

  Set_MixingParametersToDefault();

  //===============================
  // physics CP sign for aneu    //
  //===============================
  if ( NT == antineutrino ){
    for(unsigned int i = 0; i < numneu; i++){
      for(unsigned int j = i+1; j < numneu; j++){
        Set_CPPhase(i,j,-Get_CPPhase(i,j));
      }
    }
  }
  //===============================
  // init projectors             //
  //===============================

  iniProjectors();

  //===============================
  // init square mass difference //
  //===============================

  H0_array.resize(std::vector<size_t>{ne});
  for(unsigned int ie = 0; ie < ne; ie++){
    H0_array[ie] = squids::SU_vector(nsun);
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
  E_range = marray<double,1>{1};
  E_range[0] = Enu;
  Set_xrange(std::vector<double>{Enu});

  // energy is initialize
  ienergy = true;
  // state is invalited, because hamiltonian changes.
  istate = false;
}

void nuSQUIDS::init(double Emin,double Emax,unsigned int Esize, bool initialize_intereractions,double xini){
  // here the energies come in eV
  ne = Esize;
  if(elogscale){
    init(logspace(Emin,Emax,ne-1),initialize_intereractions,xini);
  }
  else{
    init(linspace(Emin,Emax,ne-1),initialize_intereractions,xini);
  }
}

void nuSQUIDS::init(marray<double,1> E_vector, bool initialize_intereractions, double xini){
  // here the energies come in natural units
  if (NT == neutrino || NT == antineutrino)
    nrhos = 1;
  else if (NT == both)
    nrhos = 2;
  else {
    throw std::runtime_error("nuSQUIDS::Error::NT = {neutrino,antineutrino,both} not : " + std::to_string(NT));
  }

  if ( numneu > SQUIDS_MAX_HILBERT_DIM )
    throw std::runtime_error("nuSQUIDS::Error::Maximum number of neutrinos exceded");
  nsun = numneu;

  ne = E_vector.size();

  //===============================
  // BEGIN                       //
  //===============================

  try {
    // initialize SQUIDS
    if (iinteraction)
      ini(ne,numneu,nrhos,nrhos,xini);
    else
      ini(ne,numneu,nrhos,0,xini);
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    throw std::runtime_error("nuSQUIDS::init : Failed while trying to initialize SQuIDS.");
  }

  try{
    SetScalarsToZero();
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    throw std::runtime_error("nuSQUIDS::init : Failed while trying to set scalar arrays to zero [SetScalarsToZero].");
  }

  Set_CoherentRhoTerms(true);
  Set_h_max(std::numeric_limits<double>::max());

  //===============================
  // initialize energy arrays    //
  //===============================
  try{
    E_range = E_vector;
    Set_xrange(std::vector<double>(E_range.begin(),E_range.end()));

    delE.resize(std::vector<size_t>{ne-1});
    for(unsigned int ei = 0; ei < ne -1; ei++)
      delE[ei] = E_range[ei+1] - E_range[ei];

    ienergy = true;
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    throw std::runtime_error("nuSQUIDS::init : Failed while constructing energy and energy difference arrays.");
  }

  //===============================
  // set parameters to default   //
  //===============================

  try{
    Set_MixingParametersToDefault();
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    throw std::runtime_error("nuSQUIDS::init : Failed while trying to set default mixing parameters [Set_MixingParametersToDetaul]");
  }

  //===============================
  // init projectors             //
  //===============================

  try{
    iniProjectors();
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    throw std::runtime_error("nuSQUIDS::init : Failed while trying initialize projectors [iniProjectors]");
  }

  //===============================
  // init square mass difference //
  //===============================

  try{
    H0_array.resize(std::vector<size_t>{ne});
    for(unsigned int ie = 0; ie < ne; ie++){
      H0_array[ie] = squids::SU_vector(nsun);
    }
    iniH0();
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    throw std::runtime_error("nuSQUIDS::init : Failed while trying initialize vaccum Hamiltonian [iniH0]");
  }

  //===============================
  // Tau properties              //
  //===============================

  taubr_lep = 0.14;
  tau_lifetime = 2.906e-13*params.sec;
  tau_mass = 1776.82*params.MeV;
  tau_reg_scale = 300.0*params.km;
  positivization_scale = 300.0*params.km;

  if(iinteraction and initialize_intereractions){
      //===============================
      // init XS and TDecay objects  //
      //===============================

      // initialize cross section object
      if ( ncs == nullptr) {
        ncs = std::make_shared<NeutrinoDISCrossSectionsFromTables>();
      } // else we assume the user has already inintialized the object if not throw error.

      try{
      // initialize tau decay spectra object
        tdc.Init(E_range[0],E_range[ne-1],ne-1);
      } catch (std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        throw std::runtime_error("nuSQUIDS::init : Failed while trying initialize TauDecaySpectra object [TauDecaySpectra::Init]");
      }
      // initialize cross section and interaction arrays
      try {
        InitializeInteractionVectors();
      } catch (std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        throw std::runtime_error("nuSQUIDS::init : Failed while trying to initialize interaction vectors [InitializeInteractionVectors]");
      }
      //===============================
      // Fill in arrays              //
      //===============================
      try {
        InitializeInteractions();
      } catch (std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        throw std::runtime_error("nuSQUIDS::init : Failed while trying to fill in interaction vectors [InitializeInteractions]");
      }
  }

  if(iinteraction){
    Set_NonCoherentRhoTerms(true);
    Set_OtherRhoTerms(true);
    Set_GammaScalarTerms(true);
    Set_OtherScalarTerms(true);
  }

  //===============================
  // END                         //
  //===============================
}

void nuSQUIDS::InitializeInteractionVectors(){

    // initialize cross section and interaction arrays
    int_struct->dNdE_NC.resize(std::vector<size_t>{nrhos,numneu,ne,ne});
    int_struct->dNdE_CC.resize(std::vector<size_t>{nrhos,numneu,ne,ne});
    int_struct->dNdE_GR.resize(std::vector<size_t>{ne,ne});
    // inverse interaction lenghts
    int_struct->invlen_NC.resize(std::vector<size_t>{nrhos,numneu,ne});
    int_struct->invlen_CC.resize(std::vector<size_t>{nrhos,numneu,ne});
    int_struct->invlen_GR.resize(std::vector<size_t>{ne});
    int_struct->invlen_INT.resize(std::vector<size_t>{nrhos,numneu,ne});
    // initialize cross section arrays
    int_struct->sigma_CC.resize(std::vector<size_t>{nrhos,numneu,ne});
    int_struct->sigma_NC.resize(std::vector<size_t>{nrhos,numneu,ne});
    int_struct->sigma_GR.resize(std::vector<size_t>{ne});
    // initialize the tau decay and interaction array
    int_struct->invlen_tau.resize(std::vector<size_t>{ne});
    int_struct->dNdE_tau_all.resize(std::vector<size_t>{ne,ne});
    int_struct->dNdE_tau_lep.resize(std::vector<size_t>{ne,ne});

}

void nuSQUIDS::PreDerive(double x){
  track->SetX(x-time_offset);
  if( basis != mass and ioscillations){
    EvolveProjectors(x);
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
  for(unsigned int rho = 0; rho < nrhos; rho++){
    for(unsigned int flv = 0; flv < numneu; flv++){
      for(unsigned int ei = 0; ei < ne; ei++){
        // will only evolve the flavor projectors
        //evol_b0_proj[rho][flv][ei] = b0_proj[flv].Evolve(h0,(x-t_ini));
        evol_b1_proj[rho][flv][ei] = b1_proj[rho][flv].Evolve(H0_array[ei],(x-Get_t_initial()));
      }
    }
  }
  return;
  /*
  if ( not use_full_hamiltonian_for_projector_evolution ){
    for(unsigned int rho = 0; rho < nrhos; rho++){
      for(unsigned int flv = 0; flv < numneu; flv++){
        for(unsigned int ei = 0; ei < ne; ei++){
          // will only evolve the flavor projectors
          //evol_b0_proj[rho][flv][ei] = b0_proj[flv].Evolve(h0,(x-t_ini));
          evol_b1_proj[rho][flv][ei] = b1_proj[rho][flv].Evolve(H0_array[ei],(x-Get_t_initial()));
        }
      }
    }
  } else {
    squids::SU_vector temp(nsun);
    for(unsigned int rho = 0; rho < nrhos; rho++){
      for(unsigned int ei = 0; ei < ne; ei++){
        temp = (H0_array[ei]+HI(ei,rho));
        temp *= (x-Get_t_initial());
        for(unsigned int flv = 0; flv < numneu; flv++){
          evol_b1_proj[rho][flv][ei] = b1_proj[rho][flv].UTransform(temp);
        }
      }
    }
  }
  */
}

squids::SU_vector nuSQUIDS::H0(double Enu, unsigned int irho) const{
  return DM2*(0.5/Enu);
}

squids::SU_vector nuSQUIDS::HI(unsigned int ie, unsigned int irho) const{
    double ye = body->ye(*track);
    double density = body->density(*track);

    double CC = params.sqrt2*params.GF*params.Na*pow(params.cm,-3)*density*ye;
    double NC;

    if (ye < 1.0e-10){
      NC = params.sqrt2*params.GF*params.Na*pow(params.cm,-3)*density;
    }
    else {
      NC = CC*(-0.5*(1.0-ye)/ye);
    }

    // construct potential in flavor basis
    //std::cout << CC << " " << NC << std::endl;
    //std::cout << irho << " " << ie << std::endl;
    //std::cout << evol_b1_proj[irho][0][ie] << std::endl;
    squids::SU_vector potential = (CC+NC)*evol_b1_proj[irho][0][ie];
    potential += (NC)*(evol_b1_proj[irho][1][ie]);
    potential += (NC)*(evol_b1_proj[irho][2][ie]);

    if ((irho == 1 and NT==both) or NT==antineutrino){
        // antineutrino matter potential flips sign
        potential *= (-1.0);
    }

    if (basis == mass){
      potential += H0_array[ie];
    }

    return potential;
    /*
    if ((irho == 0 and NT==both) or NT==neutrino){
        // neutrino potential
        return potential;
    } else if ((irho == 1 and NT==both) or NT==antineutrino){
        // antineutrino potential
        return (-1.0)*std::move(potential);
    } else{
        throw std::runtime_error("nuSQUIDS::HI : unknown particle or antiparticle");
    }
    */
}

squids::SU_vector nuSQUIDS::GammaRho(unsigned int ei,unsigned int index_rho) const{
    squids::SU_vector V(nsun);
    if (not iinteraction){
      return V;
    }

    V = evol_b1_proj[index_rho][0][ei]*(0.5*int_struct->invlen_INT[index_rho][0][ei]);
    V += evol_b1_proj[index_rho][1][ei]*(0.5*int_struct->invlen_INT[index_rho][1][ei]);
    V += evol_b1_proj[index_rho][2][ei]*(0.5*int_struct->invlen_INT[index_rho][2][ei]);

    return V;
}

squids::SU_vector nuSQUIDS::InteractionsRho(unsigned int e1,unsigned int index_rho) const{
  if(iinteraction && !ioscillations) //can use precomputed result
    return(interaction_cache[index_rho][0][e1]);

  squids::SU_vector nc_term(nsun);

  if (not iinteraction)
    return nc_term;

  // this implements the NC interactinos
  // the tau regeneration terms are implemented at the end
  // here we assume the cross section to be the same for all flavors
  squids::SU_vector projector_sum, temp;
  for(unsigned int e2 = e1 + 1; e2 < ne; e2++){
    projector_sum = evol_b1_proj[index_rho][0][e2] + evol_b1_proj[index_rho][1][e2];
    projector_sum += evol_b1_proj[index_rho][2][e2];
    temp = ACommutator(projector_sum,state[e2].rho[index_rho]);
    // NB: 0.5 compensates for factor of 2 in classical-limit {X,Y} = XY + YX = 2XY
    nc_term += temp*(0.5*int_struct->dNdE_NC[index_rho][0][e2][e1]*int_struct->invlen_NC[index_rho][0][e2])*delE[e2-1];
    // Glashow resonance
    // NB: we assume that the W- -> antineutrino branching fractions are the same for all flavors (close, but not quite)
    // to do this correctly, weight the contributions to projector_sum by the relative branching fractions
    if ((NT == both and index_rho == 1) or NT == antineutrino)
      nc_term += projector_sum
        * (evol_b1_proj[index_rho][0][e2]*state[e2].rho[index_rho])
        * (int_struct->dNdE_GR[e2][e1]*int_struct->invlen_GR[e2]*delE[e2-1]);
  }
  return nc_term;
}

double nuSQUIDS::GammaScalar(unsigned int ei, unsigned int iscalar) const{
  // we will just keep all the taus and convert them at the end
  return 0.0;
}

double nuSQUIDS::InteractionsScalar(unsigned int ei, unsigned int iscalar) const{
  if (not iinteraction)
    return 0.0;
  if(not ioscillations) //can use precomputed result
    return scalar_interaction_cache[iscalar][ei];

  double nutautoleptau = 0.0;
  for(unsigned int e2 = ei + 1; e2 < ne; e2++)
    nutautoleptau += (evol_b1_proj[iscalar][2][e2]*state[e2].rho[iscalar])*
                     (int_struct->invlen_CC[iscalar][2][e2])*(int_struct->dNdE_CC[iscalar][2][e2][ei])*delE[e2-1];
  return nutautoleptau;
}

double nuSQUIDS::GetNucleonNumber() const{
    double density = body->density(*track);
    double num_nuc = (params.gr*pow(params.cm,-3))*density*2.0/(params.proton_mass+params.neutron_mass);

    //#ifdef UpdateInteractions_DEBUG
    if(debug){
      std::cout << "============ BEGIN GetNucleonNumber ============" << std::endl;
      std::cout << "Density " << density << std::endl;
      std::cout << "Nucleon Number " << num_nuc << std::endl;
      std::cout << "============ END GetNucleonNumber ============" << std::endl;
    }
    //#endif

    if(num_nuc < 1.0e-10 ){
      num_nuc = params.Na*pow(params.cm,-3)*1.0e-10;
    }

    return num_nuc;
}

void nuSQUIDS::UpdateInteractions(){
    double num_nuc = GetNucleonNumber();
    double num_e = num_nuc*body->ye(*track);
    if(debug)
      std::cout << "============ BEGIN UpdateInteractions ============" << std::endl;
    for(unsigned int rho = 0; rho < nrhos; rho++){
      for(unsigned int flv = 0; flv < numneu; flv++){
          //#ifdef UpdateInteractions_DEBUG
          if(debug){
            std::cout << "============" << flv << "============" << std::endl;
          }
          //#endif
          for(unsigned int e1 = 0; e1 < ne; e1++){
              //#ifdef UpdateInteractions_DEBUG
              if(debug){
                  std::cout << "== CC NC Terms x = " << track->GetX()/params.km << " [km] ";
                  std::cout << "E = " << E_range[e1] << " [eV] ==" << std::endl;
                  std::cout << "CC : " << int_struct->sigma_CC[rho][flv][e1]*num_nuc << " NC : " << int_struct->sigma_NC[rho][flv][e1]*num_nuc << std::endl;
                  std::cout << "==" << std::endl;
              }
              //#endif
              int_struct->invlen_NC[rho][flv][e1] = int_struct->sigma_NC[rho][flv][e1]*num_nuc;
              int_struct->invlen_CC[rho][flv][e1] = int_struct->sigma_CC[rho][flv][e1]*num_nuc;
              int_struct->invlen_INT[rho][flv][e1] = int_struct->invlen_NC[rho][flv][e1] + int_struct->invlen_CC[rho][flv][e1];
          }
      }
    }
#ifdef WITH_GLASHOW_RESONANCE
    // Add Glashow resonance if antineutrinos are in the mix
    if (NT == both or NT == antineutrino) {
      assert(numneu > 0);
      unsigned int rho = (NT == both) ? 1 : 0;
      for(unsigned int e1 = 0; e1 < ne; e1++){
        int_struct->invlen_GR[e1] = int_struct->sigma_GR[e1]*num_e;
        int_struct->invlen_INT[rho][0][e1] += int_struct->invlen_GR[e1];
      }
    }
#endif
  
  //Without oscillations, the entries in evol_b1_proj do not depend on energy.
  //We can exploit this to precalculate the information needed by InteractionsRho
  //computing the summed projector a number of times proportional to the number
  //of energies, rather than proportinal to the square.
  if(!ioscillations){
    interaction_cache.resize(std::initializer_list<size_t>{nrhos,numneu,ne});
    squids::SU_vector projector_sum(numneu), temp(numneu);
    for(unsigned int rho = 0; rho < nrhos; rho++){
      //this will be the same for all energies, so we may as well use the first
      projector_sum = evol_b1_proj[rho][0][0] + evol_b1_proj[rho][1][0];
      projector_sum += evol_b1_proj[rho][2][0];
      
      //assume all flavors are the same, so compute for first flavor, then copy to others
      //first initialize to zero
      for(unsigned int e1=0; e1<ne; e1++){
        if(interaction_cache[rho][0][e1].Dim()!=numneu)
          interaction_cache[rho][0][e1]=squids::SU_vector(numneu);
        else
          interaction_cache[rho][0][e1].SetAllComponents(0);
      }
      //then sum contributions from higher energies cascading to lower
      for(unsigned int e2=1; e2<ne; e2++){
        temp = ACommutator(projector_sum,state[e2].rho[rho]);
        for(unsigned int e1=0; e1<e2; e1++)
          interaction_cache[rho][0][e1] += temp*(0.5*int_struct->dNdE_NC[rho][0][e2][e1]*int_struct->invlen_NC[rho][0][e2]*delE[e2-1]);
      }
      //duplicate result to other flavors
      for(unsigned int flv = 1; flv < numneu; flv++){
        for(unsigned int e1=0; e1<ne; e1++)
          interaction_cache[rho][flv][e1]=interaction_cache[rho][0][e1];
      }
#ifdef WITH_GLASHOW_RESONANCE
      // Glashow resonance
      for(unsigned int e2=1; e2<ne; e2++){
        if ((NT == both and rho == 1) or NT == antineutrino) {
          temp = projector_sum
                * (evol_b1_proj[rho][0][e2]*state[e2].rho[rho]);
          for(unsigned int e1=0; e1<e2; e1++) {
            interaction_cache[rho][0][e1] += temp*(int_struct->dNdE_GR[e2][e1]*int_struct->invlen_GR[e2]*delE[e2-1]);
          }
        }
      }
#endif
    }
    //scalar interactions for nu_tau
    scalar_interaction_cache.resize(std::initializer_list<size_t>{nrhos,ne});
    for(unsigned int rho = 0; rho < nrhos; rho++){
      //first initialize to zero
      for(unsigned int e1=0; e1<ne; e1++)
        scalar_interaction_cache[rho][e1]=0;
      for(unsigned int e2=0; e2<ne; e2++){
        double temp=(evol_b1_proj[rho][2][e2]*state[e2].rho[rho]);
        for(unsigned int e1=0; e1<e2; e1++)
          scalar_interaction_cache[rho][e1]+=temp*(int_struct->invlen_CC[rho][2][e2])*(int_struct->dNdE_CC[rho][2][e2][e1])*delE[e2-1];
      }
    }
  }
  
    if(debug)
      std::cout << "============ END UpdateInteractions ============" << std::endl;
}

void nuSQUIDS::InitializeInteractions(){

    //units
    double cm2GeV = pow(params.cm,2)*pow(params.GeV,-1);
    double cm2 = pow(params.cm,2);
    double GeVm1 = pow(params.GeV,-1);

    // load cross sections
    // initializing cross section arrays temporary array
    marray<double,4> dsignudE_CC{nrhos,numneu,ne,ne};
    marray<double,4> dsignudE_NC{nrhos,numneu,ne,ne};

    // filling cross section arrays
    std::map<unsigned int,NeutrinoCrossSections::NeutrinoType> neutype_xs_dict;
    if (NT == neutrino){
      neutype_xs_dict = (std::map<unsigned int,NeutrinoCrossSections::NeutrinoType>){{0,NeutrinoCrossSections::neutrino}};
    } else if ( NT == antineutrino ) {
      neutype_xs_dict = (std::map<unsigned int,NeutrinoCrossSections::NeutrinoType>){{0,NeutrinoCrossSections::antineutrino}};
    } else {
      // in this case NT is both
      neutype_xs_dict = (std::map<unsigned int,NeutrinoCrossSections::NeutrinoType>){{0, NeutrinoCrossSections::neutrino},{1,NeutrinoCrossSections::antineutrino}};
    }

    for(unsigned int neutype = 0; neutype < nrhos; neutype++){
      for(unsigned int flv = 0; flv < numneu; flv++){
          for(unsigned int e1 = 0; e1 < ne; e1++){
              // differential cross sections
              for(unsigned int e2 = 0; e2 < e1; e2++){
                  dsignudE_NC[neutype][flv][e1][e2] = ncs->SingleDifferentialCrossSection(E_range[e1],E_range[e2],static_cast<NeutrinoCrossSections::NeutrinoFlavor>(flv),neutype_xs_dict[neutype],NeutrinoCrossSections::NC)*cm2GeV;
                  dsignudE_CC[neutype][flv][e1][e2] = ncs->SingleDifferentialCrossSection(E_range[e1],E_range[e2],static_cast<NeutrinoCrossSections::NeutrinoFlavor>(flv),neutype_xs_dict[neutype],NeutrinoCrossSections::CC)*cm2GeV;
              }
              // total cross sections
              int_struct->sigma_CC[neutype][flv][e1] = ncs->TotalCrossSection(E_range[e1],static_cast<NeutrinoCrossSections::NeutrinoFlavor>(flv),neutype_xs_dict[neutype],NeutrinoCrossSections::CC)*cm2;
              int_struct->sigma_NC[neutype][flv][e1] = ncs->TotalCrossSection(E_range[e1],static_cast<NeutrinoCrossSections::NeutrinoFlavor>(flv),neutype_xs_dict[neutype],NeutrinoCrossSections::NC)*cm2;
          }
      }
    }

    #ifdef FixCrossSections
    // fix charge current and neutral current differential cross sections
    for(unsigned int neutype = 0; neutype < nrhos; neutype++){
      double XCC_MIN,XNC_MIN,XCC_int,XNC_int,CC_rescale,NC_rescale;
      for(unsigned int flv = 0; flv < numneu; flv++){
          XCC_MIN = int_struct->sigma_CC[neutype][flv][0];
          XNC_MIN = int_struct->sigma_CC[neutype][flv][0];
          for(unsigned int e1 = 0; e1 < ne; e1++){
              XCC_int = 0.0;
              XNC_int = 0.0;
              for(unsigned int e2 = 0; e2 < e1; e2++){
                  XCC_int += dsignudE_CC[neutype][flv][e1][e2]*delE[e2];
                  XNC_int += dsignudE_NC[neutype][flv][e1][e2]*delE[e2];
              }

              if(e1 != 0 ){
                  CC_rescale = (int_struct->sigma_CC[neutype][flv][e1] - XCC_MIN)/XCC_int;
                  NC_rescale = (int_struct->sigma_NC[neutype][flv][e1] - XNC_MIN)/XNC_int;

                  for(unsigned int e2 = 0; e2 < e1; e2++){
                      dsignudE_CC[neutype][flv][e1][e2] = dsignudE_CC[neutype][flv][e1][e2]*CC_rescale;
                      dsignudE_NC[neutype][flv][e1][e2] = dsignudE_NC[neutype][flv][e1][e2]*NC_rescale;
                  }
              }
          }
      }
    }
    #endif

    // constructing dNdE for DIS
    for(unsigned int rho = 0; rho < nrhos; rho++){
      for(unsigned int flv = 0; flv < numneu; flv++){
          for(unsigned int e1 = 0; e1 < ne; e1++){
              for(unsigned int e2 = 0; e2 < e1; e2++){
                  if (dsignudE_NC[rho][flv][e1][e2] < 1.0e-50 or (dsignudE_NC[rho][flv][e1][e2] != dsignudE_NC[rho][flv][e1][e2])){
                      int_struct->dNdE_NC[rho][flv][e1][e2] = 0.0;
                  } else {
                      int_struct->dNdE_NC[rho][flv][e1][e2] = (dsignudE_NC[rho][flv][e1][e2])/(int_struct->sigma_NC[rho][flv][e1]);
                  }
                  if (dsignudE_CC[rho][flv][e1][e2] < 1.0e-50 or (dsignudE_CC[rho][flv][e1][e2] != dsignudE_CC[rho][flv][e1][e2])){
                      int_struct->dNdE_CC[rho][flv][e1][e2] = 0.0;
                  } else {
                      int_struct->dNdE_CC[rho][flv][e1][e2] = (dsignudE_CC[rho][flv][e1][e2])/(int_struct->sigma_CC[rho][flv][e1]);
                  }
              }
          }
      }
    }
    // construct dNdE for Glashow resonance
    {
      GlashowResonanceCrossSection gr_cs;
      marray<double,2> dsignudE_GR{ne,ne};
      for(unsigned int e1 = 0; e1 < ne; e1++){
        int_struct->sigma_GR[e1] = gr_cs.TotalCrossSection(E_range[e1],NeutrinoCrossSections::electron,NeutrinoCrossSections::antineutrino,NeutrinoCrossSections::GR)*cm2;
        for(unsigned int e2 = 0; e2 < e1; e2++){
          double dsignudE = gr_cs.SingleDifferentialCrossSection(E_range[e1],E_range[e2],NeutrinoCrossSections::electron,NeutrinoCrossSections::antineutrino,NeutrinoCrossSections::GR)*cm2GeV;
          if (dsignudE < 1.0e-50 or (dsignudE != dsignudE)){
              dsignudE_GR[e1][e2] = 0.0;
          } else {
              dsignudE_GR[e1][e2] = dsignudE;
          }
        }
      }
#ifdef FixCrossSections
      // Ensure that the differential cross-sections integrate to 1
      for(unsigned int e1 = 0; e1 < ne; e1++){
        double X_int = 0.;
        for(unsigned int e2 = 0; e2 < e1; e2++){
          X_int += dsignudE_GR[e1][e2]*delE[e2];
        }
        if (e1 != 0){
          double rescale = gr_cs.GetMuonicBranchingFraction()
            * (int_struct->sigma_GR[e1] - int_struct->sigma_GR[0])/X_int;
          for(unsigned int e2 = 0; e2 < e1; e2++){
            dsignudE_GR[e1][e2] *= rescale;
          }
        }

      }
#endif
      for(unsigned int e1 = 0; e1 < ne; e1++){
        for(unsigned int e2 = 0; e2 < e1; e2++){
          int_struct->dNdE_GR[e1][e2] = dsignudE_GR[e1][e2]/int_struct->sigma_GR[e1];
        }
      }
    }

    // initialize interaction lenghts to zero
    // tau decay length array
    for(unsigned int e1 = 0; e1 < ne; e1++){
        int_struct->invlen_tau[e1] = 1.0/(tau_lifetime*E_range[e1]*tau_mass);
    }

    // load tau decay spectra

    // constructing dNdE_tau_lep/dNdE_tau_all
    for(unsigned int e1 = 0; e1 < ne; e1++){
        for(unsigned int e2 = 0; e2 < e1; e2++){
            int_struct->dNdE_tau_all[e1][e2] = tdc.dNdEnu_All(e1,e2)*GeVm1;
            int_struct->dNdE_tau_lep[e1][e2] = tdc.dNdEnu_Lep(e1,e2)*GeVm1;
        }
    }

    #ifdef FixCrossSections
    // fix tau decay spectra cross section
    double tau_all_int,tau_lep_int,tau_lep_rescale,tau_all_rescale;
    for(unsigned int e1 = 1; e1 < ne; e1++){
        tau_all_int = 0.0;
        tau_lep_int = 0.0;
        for(unsigned int e2 = 0; e2 < e1; e2++){
             tau_all_int += int_struct->dNdE_tau_all[e1][e2]*delE[e2];
             tau_lep_int += int_struct->dNdE_tau_lep[e1][e2]*delE[e2];
        }

        if( int_struct->dNdE_tau_all[e1][0]*E_range[0] < 0.25 ) {
            tau_all_rescale = (1.0 - int_struct->dNdE_tau_all[e1][0]*E_range[0])/tau_all_int;
            tau_lep_rescale = (taubr_lep - int_struct->dNdE_tau_lep[e1][0]*E_range[0])/tau_lep_int;

            for(unsigned int e2 = 0; e2 < e1; e2++){
                int_struct->dNdE_tau_all[e1][e2] = int_struct->dNdE_tau_all[e1][e2]*tau_all_rescale;
                int_struct->dNdE_tau_lep[e1][e2] = int_struct->dNdE_tau_lep[e1][e2]*tau_lep_rescale;
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
  // setting time offset
  time_offset = Get_t() - track_in->GetInitialX();
  // enforce track beginningness
  track_in->SetX(track_in->GetInitialX());
  // set track
  track = track_in;
  itrack = true;
}

void nuSQUIDS::PositivizeFlavors(){
  // advance positivity correction
  for(unsigned int rho = 0; rho < nrhos; rho++){
    for(unsigned int ie = 0; ie < ne; ie++){
      for(unsigned int flv = 0; flv < numneu; flv++){
        double quantity = EvalFlavorAtNode(flv,ie,rho);
        if( quantity < 0){
          state[ie].rho[rho] -= evol_b1_proj[rho][flv][ie]*quantity;
        }
      }
    }
  }
}

void nuSQUIDS::Set_PositivityConstrain(bool opt){
  positivization = opt;
}

void nuSQUIDS::Set_PositivityConstrainStep(double step){
  positivization_scale = step;
}

void nuSQUIDS::EvolveState(){
  // check for BODY and TRACK status
  if ( body == NULL )
    throw std::runtime_error("nuSQUIDS::Error::BODY is a NULL pointer");
  if (not ibody )
    throw std::runtime_error("nuSQUIDS::Error::Body not initialized");
  if ( track == NULL )
    throw std::runtime_error("nuSQUIDS::Error::TRACK is a NULL pointer");
  if ( not itrack )
    throw std::runtime_error("nuSQUIDS::Error::TRACK is not initialized");
  if ( not istate )
    throw std::runtime_error("nuSQUIDS::Error::Initial state not initialized");
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");

  if ( body->IsConstantDensity() and not iinteraction ){
    // when only oscillations are considered and the density is constant
    // we can ju st rotate the system to the propagation eigenstates
    use_full_hamiltonian_for_projector_evolution = true;
    // turn off squids numerics
    Set_AnyNumerics(false);
    // disable interaction basis and go to mass basis
    Set_Basis(mass);
    double evolution_time = track->GetFinalX()-track->GetInitialX();
    // go go go
    // in this case it only calls prederive and exits
    Evolve(evolution_time);

    // We will now evolve each of the states manually
    squids::SU_vector tmp1(nsun);
    squids::SU_vector tmp2(nsun);
    for(unsigned int rho = 0; rho < nrhos; rho++){
      for(unsigned int ie = 0; ie < ne; ie++){
        tmp1 = HI(ie,rho);
        //std::cout << "antes " << E_range[ie] << " " << state[ie].rho[rho] << std::endl;
        //std::cout << "hamiltonian " << tmp1 << std::endl;
        tmp2 = state[ie].rho[rho].UTransform(tmp1,gsl_complex_rect(0.,evolution_time));
        state[ie].rho[rho] = tmp2;
        //std::cout << "despues " << E_range[ie] << " " << tmp2 << std::endl;
      }
    }

    return;
  }

  if( not tauregeneration ){
    if(positivization){
      int positivization_steps = static_cast<int>((track->GetFinalX() - track->GetInitialX())/positivization_scale);
      for (int i = 0; i < positivization_steps; i++){
        Evolve(positivization_scale);
        PositivizeFlavors();
      }
      Evolve(track->GetFinalX()-positivization_scale*positivization_steps);
      PositivizeFlavors();
    } else {
      Evolve(track->GetFinalX()-track->GetInitialX());
    }
  }
  else {
    double scale;
    if(positivization)
      scale = std::min(tau_reg_scale,positivization_scale);
    else
      scale = tau_reg_scale;
    int tau_steps = static_cast<int>((track->GetFinalX() - track->GetInitialX())/scale);
    for (int i = 0; i < tau_steps; i++){
      Evolve(scale);
      if(positivization)
        PositivizeFlavors();
      ConvertTauIntoNuTau();
    }
    Evolve(track->GetFinalX()-scale*tau_steps);
    if(positivization)
      PositivizeFlavors();
    ConvertTauIntoNuTau();
  }
}

void nuSQUIDS::SetScalarsToZero(void){
  for(unsigned int rho = 0; rho < nscalars; rho++){
    for(unsigned int e1 = 0; e1 < ne; e1++){
        state[e1].scalar[rho] = 0.0;
    }
  }
}

void nuSQUIDS::ConvertTauIntoNuTau(){

  for(unsigned int e1 = 0; e1 < ne; e1++){
      double tau_neu_all  = 0.0;
      double tau_neu_lep  = 0.0;
      double tau_aneu_all = 0.0;
      double tau_aneu_lep = 0.0;

      for(unsigned int e2 = e1 +1; e2 < ne; e2++){
          //std::cout << dNdE_tau_all[e2][e1] << " " << delE[e2] << " " << state[e2].scalar[0] << std::endl;
          tau_neu_all  += int_struct->dNdE_tau_all[e2][e1]*delE[e2]*state[e2].scalar[0];
          tau_neu_lep  += int_struct->dNdE_tau_lep[e2][e1]*delE[e2]*state[e2].scalar[0];
          tau_aneu_all += int_struct->dNdE_tau_all[e2][e1]*delE[e2]*state[e2].scalar[1];
          tau_aneu_lep += int_struct->dNdE_tau_lep[e2][e1]*delE[e2]*state[e2].scalar[1];
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

void nuSQUIDS::Set_initial_state(const marray<double,1>& v, Basis basis){
  if( v.size()== 0 )
    throw std::runtime_error("nuSQUIDS::Error:Null size input array.");
  if( v.extent(0) != numneu )
    throw std::runtime_error("nuSQUIDS::Error::Initial state size not compatible with number of flavors.");
  if( not (basis == flavor || basis == mass ))
    throw std::runtime_error("nuSQUIDS::Error::BASIS can be: flavor or mass.");
  if( NT == both )
    throw std::runtime_error("nuSQUIDS::Error::Only supplied neutrino/antineutrino initial state, but set to both.");
  if( ne != 1 )
    throw std::runtime_error("nuSQUIDS::Error::nuSQUIDS initialized in multienergy mode, while state is only single energy.");
  if( !itrack or !ibody )
    throw std::runtime_error("nuSQUIDS::Error::Body and Trayectory must be specified before setting the initial state.");
  if( !ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy needs to be set before state.");

  // reset track position to initial condition
  track->SetX(track->GetInitialX());
  // we need to reinitialize the SQuIDS object to set up the new time
  ini(ne,numneu,nrhos,nscalars,track->GetInitialX());
  Set_xrange(std::vector<double>(E_range.begin(),E_range.end()));
  // since we just sync clocks, set offset to zero
  time_offset = 0;

  // initializing the projectors and hamiltonian
  SetIniFlavorProyectors();
  iniH0();

  for(unsigned int i = 0; i < ne; i++){
    for(unsigned int r = 0; r < nrhos; r++){
      if (basis == flavor){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(unsigned int j = 0; j < v.extent(0); j++)
        {
          state[i].rho[r] += v[j]*b1_proj[r][j];
        }
      }
      else if (basis == mass){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(int j = 0; j < v.extent(0); j++)
        {
          state[i].rho[r] += v[j]*b0_proj[j];
        }
      }
    }
  }
  if(nscalars)
    SetScalarsToZero();

  istate = true;
};

void nuSQUIDS::Set_initial_state(const marray<double,2>& v, Basis basis){
  if( v.size() == 0 )
    throw std::runtime_error("nuSQUIDS::Error:Null size input array.");
  if( v.extent(0) != ne )
    throw std::runtime_error("nuSQUIDS::Error:Input vector with wrong dimensions.(ne:"+std::to_string(v.extent(0))+"!="+std::to_string(ne)+")");
  if( v.extent(1) != numneu)
    throw std::runtime_error("nuSQUIDS::Error:Input vector with wrong dimensions.(numneu:"+std::to_string(v.extent(1))+"!="+std::to_string(numneu)+")");
  if( not (basis == flavor || basis == mass ))
    throw std::runtime_error("nuSQUIDS::Error::BASIS can be : flavor or mass.");
  if( NT == both )
    throw std::runtime_error("nuSQUIDS::Error::Only supplied neutrino/antineutrino initial state, but set to both.");
  if( !itrack or !ibody )
    throw std::runtime_error("nuSQUIDS::Error::Body and Trayectory must be specified before setting the initial state.");
  if( !ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy needs to be set before state.");

  // reset track position to initial condition
  track->SetX(track->GetInitialX());
  // we need to reinitialize the SQuIDS object to set up the new time
  ini(ne,numneu,nrhos,nscalars,track->GetInitialX());
  Set_xrange(std::vector<double>(E_range.begin(),E_range.end()));
  // since we just sync clocks, set offset to zero
  time_offset = 0;

  // initializing the projectors and hamiltonian
  SetIniFlavorProyectors();
  iniH0();

  for(unsigned int i = 0; i < ne; i++){
    for(unsigned int r = 0; r < nrhos; r++){
      if (basis == flavor){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(unsigned int j = 0; j < numneu; j++)
        {
          state[i].rho[r] += v[i][j]*b1_proj[r][j];
        }
      }
      else if (basis == mass){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(unsigned int j = 0; j < numneu; j++){
          state[i].rho[r] += v[i][j]*b0_proj[j];
        }
      }
    }
  }
  if(nscalars)
    SetScalarsToZero();

  istate = true;
}

void nuSQUIDS::Set_initial_state(const marray<double,3>& v, Basis basis){
  if( v.size() == 0 )
    throw std::runtime_error("nuSQUIDS::Error:Null size input array.");
  if( v.extent(0) != ne )
    throw std::runtime_error("nuSQUIDS::Error:Input vector with wrong dimensions.(ne:"+std::to_string(v.extent(0))+"!="+std::to_string(ne)+")");
  if( v.extent(1) != nrhos)
    throw std::runtime_error("nuSQUIDS::Error:Input vector with wrong dimensions.(nrhos:"+std::to_string(v.extent(1))+"!="+std::to_string(nrhos)+")");
  if( v.extent(2) != numneu)
    throw std::runtime_error("nuSQUIDS::Error:Input vector with wrong dimensions.(numneu"+std::to_string(v.extent(2))+"!="+std::to_string(numneu)+")");
  if( not (basis == flavor || basis == mass ))
    throw std::runtime_error("nuSQUIDS::Error::BASIS can be : flavor or mass.");
  if( NT != both )
    throw std::runtime_error("nuSQUIDS::Error::Supplied neutrino and antineutrino initial state, but not set to both.");
  if( !itrack or !ibody )
    throw std::runtime_error("nuSQUIDS::Error::Body and Trayectory must be specified before setting the initial state.");
  if( !ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy needs to be set before state.");

  // reset track position to initial condition
  track->SetX(track->GetInitialX());
  // we need to reinitialize the SQuIDS object to set up the new time
  ini(ne,numneu,nrhos,nscalars,track->GetInitialX());
  // since we just sync clocks, set offset to zero
  time_offset = 0;

  Set_xrange(std::vector<double>(E_range.begin(),E_range.end()));

  // initializing the projectors and hamiltonian
  SetIniFlavorProyectors();
  if (!ienergy)
    std::cout << "energy not set" << std::endl;
  iniH0();

  for(unsigned int i = 0; i < ne; i++){
    for(unsigned int r = 0; r < nrhos; r++){
      if (basis == flavor){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(unsigned int j = 0; j < numneu; j++)
        {
          state[i].rho[r] += v[i][r][j]*b1_proj[r][j];
        }
      }
      else if (basis == mass){
        state[i].rho[r] = 0.0*b0_proj[0];
        for(unsigned int j = 0; j < numneu; j++){
          state[i].rho[r] += v[i][r][j]*b0_proj[j];
        }
      }
    }
  }
  if(nscalars)
    SetScalarsToZero();
  istate = true;
}

marray<double,1> nuSQUIDS::GetERange() const{
  return E_range;
}

size_t nuSQUIDS::GetNumE() const{
  return ne;
}

double nuSQUIDS::EvalMass(unsigned int flv,double EE, unsigned int rho) const{
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != both )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  if ( basis == mass )
    throw std::runtime_error("nuSQUIDS::Error::Use EvalMassAtNode. Interpolation is not recommended on this basis.");
  return GetExpectationValueD(b0_proj[flv], rho, EE);
}

double nuSQUIDS::EvalFlavor(unsigned int flv,double EE,unsigned int rho) const{
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != both )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  if ( basis == mass )
    throw std::runtime_error("nuSQUIDS::Error::Use EvalMassAtNode. Interpolation is not recommended on this basis.");
  return GetExpectationValueD(b1_proj[rho][flv], rho, EE);
}

double nuSQUIDS::EvalMassAtNode(unsigned int flv, unsigned int ei, unsigned int rho) const{
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != both )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  if(basis == mass)
    return b0_proj[flv]*state[ei].rho[rho];
  return GetExpectationValue(b0_proj[flv], rho, ei);
}

double nuSQUIDS::EvalFlavorAtNode(unsigned int flv, unsigned int ei, unsigned int rho) const{
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != both )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  if(basis == mass)
    return b1_proj[rho][flv]*state[ei].rho[rho];
  return GetExpectationValue(b1_proj[rho][flv], rho, ei);
}

double nuSQUIDS::EvalMass(unsigned int flv) const{
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

double nuSQUIDS::EvalFlavor(unsigned int flv) const{
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
  if(use_full_hamiltonian_for_projector_evolution)
    return evol_b1_proj[0][flv][0]*state[0].rho[0];
  return GetExpectationValue(b1_proj[0][flv], 0, 0);
}

void nuSQUIDS::iniH0(){
  DM2 = squids::SU_vector(nsun);
  for(unsigned int i = 1; i < nsun; i++){
      DM2 += (b0_proj[i])*params.GetEnergyDifference(i);
  }

  if(ienergy){
    for(unsigned int ei = 0; ei < ne; ei++){
      // here we assume that H0 is the same for neutrinos and antineutrinos
      // else we will need two H0_array.
      H0_array[ei] = H0(E_range[ei],0);
    }
  }
}

void nuSQUIDS::AntineutrinoCPFix(unsigned int rho){
  if(NT == antineutrino or (NT == both and rho == 1)){
    for(unsigned int i = 0; i < numneu; i++){
      for(unsigned int j = i+1; j < numneu; j++){
        Set_CPPhase(i,j,-Get_CPPhase(i,j));
      }
    }
  }
}

void nuSQUIDS::iniProjectors(){

  b0_proj.resize(std::vector<size_t>{numneu});
  for(unsigned int flv = 0; flv < numneu; flv++){
    b0_proj[flv] = squids::SU_vector::Projector(nsun,flv);
  }

  b1_proj.resize(std::vector<size_t>{nrhos,numneu});
  for(unsigned int rho = 0; rho < nrhos; rho++){
    for(unsigned int flv = 0; flv < numneu; flv++){
      b1_proj[rho][flv] = squids::SU_vector::Projector(nsun,flv);

      AntineutrinoCPFix(rho);
      b1_proj[rho][flv].RotateToB1(params);
      AntineutrinoCPFix(rho);
    }
  }

  evol_b0_proj.resize(std::vector<size_t>{nrhos,numneu,ne});
  evol_b1_proj.resize(std::vector<size_t>{nrhos,numneu,ne});
  for(unsigned int rho = 0; rho < nrhos; rho++){
    for(unsigned int flv = 0; flv < numneu; flv++){
      for(unsigned int e1 = 0; e1 < ne; e1++){
        evol_b0_proj[rho][flv][e1] = squids::SU_vector::Projector(nsun,flv);
        evol_b1_proj[rho][flv][e1] = squids::SU_vector::Projector(nsun,flv);

        AntineutrinoCPFix(rho);
        evol_b1_proj[rho][flv][e1].RotateToB1(params);
        AntineutrinoCPFix(rho);
      }
    }
  }

}

void nuSQUIDS::SetIniFlavorProyectors(){
  for(unsigned int rho = 0; rho < nrhos; rho++){
    for(unsigned int flv = 0; flv < numneu; flv++){
      for(unsigned int e1 = 0; e1 < ne; e1++){
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

const squids::SU_vector& nuSQUIDS::GetState(unsigned int ie, unsigned int rho) const{
  return state[ie].rho[rho];
}

const squids::SU_vector& nuSQUIDS::GetFlavorProj(unsigned int flv,unsigned int rho) const{
  return b1_proj[rho][flv];
}

const squids::SU_vector& nuSQUIDS::GetMassProj(unsigned int flv,unsigned int rho) const{
  return b0_proj[flv];
}

squids::SU_vector nuSQUIDS::GetHamiltonian(unsigned int ei, unsigned int rho){
  if (!ienergy)
    throw std::runtime_error("nuSQUIDS::Error::Energy not initialized");
  PreDerive(Get_t());
  return H0(E_range[ei],rho)+HI(ei,rho,Get_t());
}

void nuSQUIDS::WriteStateHDF5(std::string str,std::string grp,bool save_cross_section, std::string cross_section_grp_loc) const{
  if ( body == NULL )
    throw std::runtime_error("nuSQUIDS::Error::BODY is a NULL pointer");
  if (not ibody )
    throw std::runtime_error("nuSQUIDS::Error::Body not initialized");
  if ( track == NULL )
    throw std::runtime_error("nuSQUIDS::Error::TRACK is a NULL pointer");
  if ( not itrack )
    throw std::runtime_error("nuSQUIDS::Error::TRACK is not initialized");
  if ( not istate )
    throw std::runtime_error("nuSQUIDS::Error::Initial state not initialized");
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");

  if (!iinteraction)
    save_cross_section = iinteraction;

  // this lines supress HDF5 error messages
  H5Eset_auto (H5E_DEFAULT,NULL, NULL);

  hid_t file_id,group_id,root_id;
  hid_t dset_id;
  // create HDF5 file
  //std::cout << "writing to hdf5 file" << std::endl;
  // H5F_ACC_TRUNC : overwrittes file
  // H5F_ACC_EXCL  : files if file exists
  file_id = H5Fopen(str.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0 ) {// file already exists
    file_id = H5Fcreate(str.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
        throw std::runtime_error("nuSQUIDS::Error::Cannot create file at " + str + ".");
  }
  root_id = H5Gopen(file_id, "/",H5P_DEFAULT);
  if (strncmp(grp.c_str(),"/",1)!=0){
    std::cout << "nuSQUIDS::WriteStateHDF5::Warning::group location name did not start with '/'. '/' was prepend to " + grp << std::endl;
    grp = "/"+grp;
  }

  if ( grp != "/" )
    group_id = H5Gcreate(root_id, grp.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  else
    group_id = root_id;

  // write the energy range
  hsize_t Edims[1]={E_range.extent(0)};
  dset_id = H5LTmake_dataset(group_id,"energies",1,Edims,H5T_NATIVE_DOUBLE,E_range.get_data());
  H5LTset_attribute_string(group_id, "energies", "elogscale", (elogscale) ? "True":"False");

  // write mixing parameters
  hsize_t dim[1]{1};
  H5LTmake_dataset(group_id,"basic",1,dim,H5T_NATIVE_DOUBLE,0);
  H5LTmake_dataset(group_id,"mixingangles",1,dim,H5T_NATIVE_DOUBLE,0);
  H5LTmake_dataset(group_id,"CPphases",1,dim,H5T_NATIVE_DOUBLE,0);
  H5LTmake_dataset(group_id,"massdifferences",1,dim,H5T_NATIVE_DOUBLE,0);

  H5LTset_attribute_int(group_id, "basic","numneu",(const int*)&numneu, 1);
  int auxint = static_cast<int>(NT);
  H5LTset_attribute_int(group_id, "basic","NT",&auxint,1);
  H5LTset_attribute_string(group_id, "basic", "interactions", (iinteraction) ? "True":"False");
  double auxt = Get_t();
  H5LTset_attribute_double(group_id, "basic", "squids_time", &auxt,1);
  double auxt_ini = Get_t_initial();
  H5LTset_attribute_double(group_id, "basic", "squids_time_initial", &auxt_ini,1);

  // version numbers
  H5LTset_attribute_string(group_id, "basic", "squids_version", SQUIDS_VERSION_STR);
  unsigned int squids_version = SQUIDS_VERSION;
  H5LTset_attribute_uint(group_id, "basic", "squids_version_number", &squids_version,1);

  H5LTset_attribute_string(group_id, "basic", "nusquids_version", NUSQUIDS_VERSION_STR);
  unsigned int nusquids_version = NUSQUIDS_VERSION;
  H5LTset_attribute_uint(group_id, "basic", "nusquids_version_number", &nusquids_version,1);

  // set mixing angles
  for( unsigned int i = 0; i < numneu; i++ ){
    for( unsigned int j = i+1; j < numneu; j++ ){
      std::string th_label = "th"+std::to_string(i+1)+std::to_string(j+1);
      double th_value = params.GetMixingAngle(i,j);
      H5LTset_attribute_double(group_id, "mixingangles",th_label.c_str(),&th_value, 1);

      std::string delta_label = "delta"+std::to_string(i+1)+std::to_string(j+1);
      double delta_value = params.GetPhase(i,j);
      H5LTset_attribute_double(group_id, "CPphases",delta_label.c_str(),&delta_value, 1);
    }
  }

  for ( unsigned int i = 1; i < numneu; i++ ){
    std::string dm2_label = "dm"+std::to_string(i+1)+"1sq";
    double dm2_value = params.GetEnergyDifference(i);
    H5LTset_attribute_double(group_id, "massdifferences",dm2_label.c_str(),&dm2_value, 1);
  }

  //writing state
  const unsigned int numneusq = numneu*numneu;
  hsize_t statedim[2] {E_range.size(),(hsize_t)numneu*numneu};
  std::vector<double> neustate(numneusq*ne), aneustate(numneusq*ne);

  for(unsigned int ie = 0; ie < ne; ie++){
    for(unsigned int i = 0; i < numneu*numneu; i ++){
      if (NT == both){
        neustate[ie*numneusq + i] = state[ie].rho[0][i];
        aneustate[ie*numneusq + i] = state[ie].rho[1][i];
      }
      else if (NT == neutrino){
        neustate[ie*numneusq + i] = state[ie].rho[0][i];
        aneustate[ie*numneusq + i] = 0.0;
      }
      else if (NT == antineutrino){
        neustate[ie*numneusq + i] = 0.0;
        aneustate[ie*numneusq + i] = state[ie].rho[0][i];
      }
    }
  }

  dset_id = H5LTmake_dataset(group_id,"neustate",2,statedim,H5T_NATIVE_DOUBLE,static_cast<const void*>(neustate.data()));
  dset_id = H5LTmake_dataset(group_id,"aneustate",2,statedim,H5T_NATIVE_DOUBLE,static_cast<const void*>(aneustate.data()));

/*
  // writing state flavor and mass composition
  hsize_t pdim[2] {E_range.size(), static_cast<hsize_t>(numneu)};
  if ( NT == both )
    pdim[1] *= 2;
  std::vector<double> flavor,mass;

  for(unsigned int ie = 0; ie < ne; ie++){
    // neutrino
    if (NT == both or NT == neutrino){
      for(unsigned int i = 0; i < numneu; i++){
          flavor.push_back(EvalFlavorAtNode(i,ie,0));
          mass.push_back(EvalMassAtNode(i,ie,0));
        }
    }
      // antineutrino
    if (NT == both or NT == antineutrino){
      for(unsigned int i = 0; i < numneu; i++){
          flavor.push_back(EvalFlavorAtNode(i,ie,0));
          mass.push_back(EvalMassAtNode(i,ie,0));
      }
    }
  }

  dset_id = H5LTmake_dataset(group_id,"flavorcomp",2,pdim,H5T_NATIVE_DOUBLE,static_cast<const void*>(flavor.data()));
  idset_id = H5LTmake_dataset(group_id,"masscomp",2,pdim,H5T_NATIVE_DOUBLE,static_cast<const void*>(mass.data()));
*/

  // writing body and track information
  hsize_t trackparamdim[1] {track->GetTrackParams().size()};
  if ( trackparamdim[0] == 0 ) {
    H5LTmake_dataset(group_id,"track",1,dim,H5T_NATIVE_DOUBLE,0);
  } else {
    H5LTmake_dataset(group_id,"track",1,trackparamdim,H5T_NATIVE_DOUBLE,static_cast<const void*>(track->GetTrackParams().data()));
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
    H5LTmake_dataset(group_id,"body",1,bodyparamdim,H5T_NATIVE_DOUBLE,static_cast<const void*>(body->GetBodyParams().data()));
  }
  H5LTset_attribute_string(group_id, "body", "NAME", body->GetName().c_str());
  unsigned int bid = body->GetId();
  H5LTset_attribute_uint(group_id, "body", "ID", &bid,1);

  // writing cross section information
  hid_t xs_group_id;
  if ( cross_section_grp_loc == ""){
    xs_group_id = H5Gcreate(group_id, "crosssections", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  } else {
    xs_group_id = H5Gcreate(root_id, cross_section_grp_loc.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  if (iinteraction and save_cross_section) {
    // sigma_CC and sigma_NC
    hsize_t XSdim[3] {static_cast<hsize_t>(nrhos),
                      static_cast<hsize_t>(numneu),
                      static_cast<hsize_t>(ne)};
    std::vector<double> xsCC(nrhos*numneu*ne),xsNC(nrhos*numneu*ne),xsGR(ne);
    for ( unsigned int rho = 0; rho < nrhos; rho ++){
      for ( unsigned int flv = 0; flv < numneu; flv ++){
          for ( unsigned int ie = 0; ie < ne; ie ++){
            xsCC[rho*(numneu*ne) +  flv*ne + ie] = int_struct->sigma_CC[rho][flv][ie];
            xsNC[rho*(numneu*ne) +  flv*ne + ie] = int_struct->sigma_NC[rho][flv][ie];
          }
      }
    }
    std::copy(int_struct->sigma_GR.begin(), int_struct->sigma_GR.end(), xsGR.begin());
    dset_id = H5LTmake_dataset(xs_group_id,"sigmacc",3,XSdim,H5T_NATIVE_DOUBLE,static_cast<const void*>(xsCC.data()));
    dset_id = H5LTmake_dataset(xs_group_id,"sigmanc",3,XSdim,H5T_NATIVE_DOUBLE,static_cast<const void*>(xsNC.data()));
    dset_id = H5LTmake_dataset(xs_group_id,"sigmagr",1,&XSdim[2],H5T_NATIVE_DOUBLE,static_cast<const void*>(xsGR.data()));

    // dNdE_CC and dNdE_NC
    hsize_t dXSdim[4] {static_cast<hsize_t>(nrhos),
                       static_cast<hsize_t>(numneu),
                       static_cast<hsize_t>(ne),
                       static_cast<hsize_t>(ne)};
    std::vector<double> dxsCC(nrhos*numneu*ne*ne),dxsNC(nrhos*numneu*ne*ne),dxsGR(ne*ne);

    for(unsigned int rho = 0; rho < nrhos; rho++){
      for(unsigned int flv = 0; flv < numneu; flv++){
          for(unsigned int e1 = 0; e1 < ne; e1++){
              for(unsigned int e2 = 0; e2 < ne; e2++){
                if (e2 < e1) {
                  dxsCC[rho*(numneu*ne*ne) +  flv*ne*ne + e1*ne + e2] = int_struct->dNdE_CC[rho][flv][e1][e2];
                  dxsNC[rho*(numneu*ne*ne) +  flv*ne*ne + e1*ne + e2] = int_struct->dNdE_NC[rho][flv][e1][e2];
                } else {
                  dxsCC[rho*(numneu*ne*ne) +  flv*ne*ne + e1*ne + e2] = 0.0;
                  dxsNC[rho*(numneu*ne*ne) +  flv*ne*ne + e1*ne + e2] = 0.0;
                }
              }
          }
      }
    }
    std::copy(int_struct->dNdE_GR.begin(), int_struct->dNdE_GR.end(), dxsGR.begin());
    dset_id = H5LTmake_dataset(xs_group_id,"dNdEcc",4,dXSdim,H5T_NATIVE_DOUBLE,static_cast<const void*>(dxsCC.data()));
    dset_id = H5LTmake_dataset(xs_group_id,"dNdEnc",4,dXSdim,H5T_NATIVE_DOUBLE,static_cast<const void*>(dxsNC.data()));
    dset_id = H5LTmake_dataset(xs_group_id,"dNdEgr",2,&dXSdim[2],H5T_NATIVE_DOUBLE,static_cast<const void*>(dxsGR.data()));

    // invlen_tau
    hsize_t iltdim[1] {static_cast<hsize_t>(ne)};
    dset_id = H5LTmake_dataset(xs_group_id,"invlentau",1,iltdim,H5T_NATIVE_DOUBLE,static_cast<const void*>(int_struct->invlen_tau.get_data()));

    // dNdE_tau_all,dNdE_tau_lep
    hsize_t dNdEtaudim[2] {static_cast<hsize_t>(ne),
                           static_cast<hsize_t>(ne)};
    std::vector<double> dNdEtauall(ne*ne),dNdEtaulep(ne*ne);
    for(unsigned int e1 = 0; e1 < ne; e1++){
        for(unsigned int e2 = 0; e2 < ne; e2++){
          if ( e2 < e1 ) {
            dNdEtauall[e1*ne + e2] = int_struct->dNdE_tau_all[e1][e2];
            dNdEtaulep[e1*ne + e2] = int_struct->dNdE_tau_lep[e1][e2];
          } else  {
            dNdEtauall[e1*ne + e2] = 0.0;
            dNdEtaulep[e1*ne + e2] = 0.0;
          }
        }
    }

    dset_id = H5LTmake_dataset(xs_group_id,"dNdEtauall",2,dNdEtaudim,H5T_NATIVE_DOUBLE,static_cast<void*>(dNdEtauall.data()));
    dset_id = H5LTmake_dataset(xs_group_id,"dNdEtaulep",2,dNdEtaudim,H5T_NATIVE_DOUBLE,static_cast<void*>(dNdEtaulep.data()));
  }

  // close cross section group
  H5Gclose(xs_group_id);

  // write user parameters
  hid_t user_parameters_id = H5Gcreate(group_id, "user_parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // give control to the user and temporary restore HDF5 error messages
  //H5Eset_auto (H5E_DEFAULT,(H5E_auto_t) H5Eprint,stderr);
  AddToWriteHDF5(user_parameters_id);
  //H5Eset_auto (H5E_DEFAULT,NULL, NULL);

  // close root group
  H5Gclose ( root_id );
  if ( root_id != group_id )
    H5Gclose ( group_id );
  // close HDF5 file
  H5Fclose (file_id);
}

void nuSQUIDS::AddToWriteHDF5(hid_t hdf5_loc_id) const {

}

void nuSQUIDS::AddToReadHDF5(hid_t hdf5_loc_id){

}

void nuSQUIDS::ReadStateHDF5(std::string str,std::string grp,std::shared_ptr<InteractionStructure> iis){
  hid_t file_id,group_id,root_id;
  // open HDF5 file
  //std::cout << "reading from hdf5 file" << std::endl;
  file_id = H5Fopen(str.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
      throw std::runtime_error("nuSQUIDS::Error::file not found : " + str + ".");
  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  if (strncmp(grp.c_str(),"/",1)!=0){
    std::cout << "nuSQUIDS::ReadStateHDF5::Warning::group location name did not start with '/'. '/' was prepend to " + grp << std::endl;
    grp = "/"+grp;
  }
  group_id = H5Gopen(root_id, grp.c_str(), H5P_DEFAULT);
  if ( group_id < 0 )
      throw std::runtime_error("nuSQUIDS::Error::Group '" + grp + "' does not exist in HDF5.");

  // read number of neutrinos
  H5LTget_attribute_uint(group_id, "basic", "numneu", static_cast<unsigned int*>(&numneu));
  // neutrino/antineutrino/both
  int auxint;
  char auxchar[20];
  H5LTget_attribute_int(group_id, "basic", "NT", &auxint);
  NT = static_cast<NeutrinoType>(auxint);
  // interactions
  H5LTget_attribute_string(group_id,"basic","interactions", auxchar);
  std::string aux = auxchar;
  if ( aux == "True")
    iinteraction = true;
  else
    iinteraction = false;

  double squids_time;
  H5LTget_attribute_double(group_id, "basic", "squids_time", &squids_time);

  double squids_time_initial;
  H5LTget_attribute_double(group_id, "basic", "squids_time_initial", &squids_time_initial);

  // check version numbers
  unsigned int squids_version;
  H5LTget_attribute_uint(group_id, "basic", "squids_version_number", &squids_version);

  if ( squids_version > SQUIDS_VERSION )
    throw std::runtime_error("nuSQUIDS::ReadStateHDF5::Error: File was written using SQuIDS version " +
        std::to_string(squids_version) + " current version is " + std::to_string(SQUIDS_VERSION));

  unsigned int nusquids_version;
  H5LTget_attribute_uint(group_id, "basic", "nusquids_version_number", &nusquids_version);
  if ( nusquids_version > NUSQUIDS_VERSION )
    throw std::runtime_error("nuSQUIDS::ReadStateHDF5::Error: File was written using nuSQuIDS version " +
        std::to_string(nusquids_version) + " current version is " + std::to_string(NUSQUIDS_VERSION));

  // read and set mixing parameters
  for( unsigned int i = 0; i < numneu; i++ ){
    for( unsigned int j = i+1; j < numneu; j++ ){
      double th_value;
      std::string th_label = "th"+std::to_string(i+1)+std::to_string(j+1);
      H5LTget_attribute_double(group_id,"mixingangles", th_label.c_str(), &th_value);
      Set_MixingAngle(i,j,th_value);

      double delta_value;
      std::string delta_label = "delta"+std::to_string(i+1)+std::to_string(j+1);
      H5LTget_attribute_double(group_id,"CPphases", delta_label.c_str(), &delta_value);
      Set_CPPhase(i,j,delta_value);
    }
  }

  for( unsigned int i = 1; i < numneu; i++ ){
    double dm2_value;
    std::string dm2_label = "dm"+std::to_string(i+1)+"1sq";
    H5LTget_attribute_double(group_id,"massdifferences", dm2_label.c_str(), &dm2_value);
    Set_SquareMassDifference(i, dm2_value);
  }

  // reading energy
  hsize_t dims[2];
  H5LTget_dataset_info(group_id, "energies", dims, NULL, NULL);

  ne = static_cast<unsigned int>(dims[0]);
  std::unique_ptr<double[]> energy_data(new double[ne]);
  H5LTread_dataset_double(group_id, "energies", energy_data.get());
  E_range = marray<double,1>{ne};
  for (unsigned int ie = 0; ie < ne; ie++)
    E_range[ie] = energy_data[ie];

  H5LTget_attribute_string(group_id,"energies","elogscale", auxchar);
  aux = auxchar;
  if ( aux == "True")
    elogscale = true;
  else
    elogscale = false;

  // reading body and track
  unsigned int body_id;
  hsize_t dimbody[1];
  H5LTget_attribute_uint(group_id,"body","ID",&body_id);

  H5LTget_dataset_info(group_id,"body", dimbody,NULL,NULL);
  double body_params[dimbody[0]];
  H5LTread_dataset_double(group_id,"body", body_params);

  hsize_t dimtrack[1];
  H5LTget_dataset_info(group_id,"track", dimtrack ,NULL,NULL);
  std::unique_ptr<double[]> track_params(new double[dimtrack[0]]);
  H5LTread_dataset_double(group_id,"track", track_params.get());

  double x_current;
  H5LTget_attribute_double(group_id,"track","X",&x_current);

  // setting body and track
  SetBodyTrack(body_id,dimbody[0],body_params,dimtrack[0],track_params.get());

  // set trayectory to current time
  track->SetX(x_current);

  // initializing nuSQUIDS
  if (ne == 1){
    if(not inusquids)
      init(squids_time_initial);
    Set_E(energy_data[0]);
  }
  else {
    init(E_range,false,squids_time_initial);
  }
  // reset current squids time
  Set_t(squids_time);
  // set time offset
  time_offset = squids_time - track->GetX();

  // evolve projectors to current time
  EvolveProjectors(squids_time);
  // reading state
  H5LTget_dataset_info(group_id,"neustate", dims,NULL,NULL);
  std::unique_ptr<double[]> neudata(new double[dims[0]*dims[1]]);
  H5LTread_dataset_double(group_id,"neustate", neudata.get());

  H5LTget_dataset_info(group_id,"aneustate", dims,NULL,NULL);
  std::unique_ptr<double[]> aneudata(new double[dims[0]*dims[1]]);
  H5LTread_dataset_double(group_id,"aneustate", aneudata.get());

  for(unsigned int ie = 0; ie < dims[0]; ie++){
    for (unsigned int j = 0; j < dims[1]; j ++){
      if (NT == neutrino)
        state[ie].rho[0][j] = neudata[ie*dims[1]+j];
      else if ( NT == antineutrino)
        state[ie].rho[0][j] = aneudata[ie*dims[1]+j];
      else if ( NT == both ){
        state[ie].rho[0][j] = neudata[ie*dims[1]+j];
        state[ie].rho[1][j] = aneudata[ie*dims[1]+j];
      }
    }
  }

  if(iinteraction) {
    if ( iis == nullptr ){
      throw std::runtime_error("nuSQUIDS::ReadStateHDF5::No interaction structure provided.");
    } else {
      int_struct = iis;
    }
  }

  // read from user parameters
  hid_t user_parameters_id = H5Gopen(group_id, "user_parameters", H5P_DEFAULT);
  //H5Eset_auto (H5E_DEFAULT,(H5E_auto_t) H5Eprint,stderr);
  AddToReadHDF5(user_parameters_id);
  //H5Eset_auto (H5E_DEFAULT,NULL, NULL);
  H5Gclose(user_parameters_id);

  // close HDF5 file
  H5Gclose ( group_id );
  // close root and file
  H5Gclose ( root_id );
  H5Fclose (file_id);

  // we assume that this was created with the writer and got to this point!
  istate = true;
  ienergy = true;
  itrack = true;
  ibody = true;
  // initialize H0
  iniH0();
}

void nuSQUIDS::ReadStateHDF5(std::string str,std::string grp,std::string cross_section_grp_loc){
  hid_t file_id,group_id,root_id;
  // open HDF5 file
  //std::cout << "reading from hdf5 file" << std::endl;
  file_id = H5Fopen(str.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
      throw std::runtime_error("nuSQUIDS::Error::file not found : " + str + ".");
  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  if (strncmp(grp.c_str(),"/",1)!=0){
    std::cout << "nuSQUIDS::ReadStateHDF5::Warning::group location name did not start with '/'. '/' was prepend to " + grp << std::endl;
    grp = "/"+grp;
  }
  group_id = H5Gopen(root_id, grp.c_str(), H5P_DEFAULT);
  if ( group_id < 0 )
      throw std::runtime_error("nuSQUIDS::Error::Group '" + grp + "' does not exist in HDF5.");

  // read number of neutrinos
  H5LTget_attribute_uint(group_id, "basic", "numneu", static_cast<unsigned int*>(&numneu));
  // neutrino/antineutrino/both
  int auxint;
  char auxchar[20];
  H5LTget_attribute_int(group_id, "basic", "NT", &auxint);
  NT = static_cast<NeutrinoType>(auxint);
  // interactions
  H5LTget_attribute_string(group_id,"basic","interactions", auxchar);
  std::string aux = auxchar;
  if ( aux == "True")
    iinteraction = true;
  else
    iinteraction = false;

  double squids_time;
  H5LTget_attribute_double(group_id, "basic", "squids_time", &squids_time);

  double squids_time_initial;
  H5LTget_attribute_double(group_id, "basic", "squids_time_initial", &squids_time_initial);

  // check version numbers
  unsigned int squids_version;
  H5LTget_attribute_uint(group_id, "basic", "squids_version_number", &squids_version);

  if ( squids_version > SQUIDS_VERSION )
    throw std::runtime_error("nuSQUIDS::ReadStateHDF5::Error: File was written using SQuIDS version " +
        std::to_string(squids_version) + " current version is " + std::to_string(SQUIDS_VERSION));

  unsigned int nusquids_version;
  H5LTget_attribute_uint(group_id, "basic", "nusquids_version_number", &nusquids_version);
  if ( nusquids_version > NUSQUIDS_VERSION )
    throw std::runtime_error("nuSQUIDS::ReadStateHDF5::Error: File was written using nuSQuIDS version " +
        std::to_string(nusquids_version) + " current version is " + std::to_string(NUSQUIDS_VERSION));

  // read and set mixing parameters
  for( unsigned int i = 0; i < numneu; i++ ){
    for( unsigned int j = i+1; j < numneu; j++ ){
      double th_value;
      std::string th_label = "th"+std::to_string(i+1)+std::to_string(j+1);
      H5LTget_attribute_double(group_id,"mixingangles", th_label.c_str(), &th_value);
      Set_MixingAngle(i,j,th_value);

      double delta_value;
      std::string delta_label = "delta"+std::to_string(i+1)+std::to_string(j+1);
      H5LTget_attribute_double(group_id,"CPphases", delta_label.c_str(), &delta_value);
      Set_CPPhase(i,j,delta_value);
    }
  }

  for( unsigned int i = 1; i < numneu; i++ ){
    double dm2_value;
    std::string dm2_label = "dm"+std::to_string(i+1)+"1sq";
    H5LTget_attribute_double(group_id,"massdifferences", dm2_label.c_str(), &dm2_value);
    Set_SquareMassDifference(i, dm2_value);
  }

  // reading energy
  hsize_t dims[2];
  H5LTget_dataset_info(group_id, "energies", dims, NULL, NULL);

  ne = static_cast<unsigned int>(dims[0]);
  std::unique_ptr<double[]> energy_data(new double[ne]);
  H5LTread_dataset_double(group_id, "energies", energy_data.get());
  E_range = marray<double,1>{ne};
  for (unsigned int ie = 0; ie < ne; ie++)
    E_range[ie] = energy_data[ie];

  H5LTget_attribute_string(group_id,"energies","elogscale", auxchar);
  aux = auxchar;
  if ( aux == "True")
    elogscale = true;
  else
    elogscale = false;

  // reading body and track
  unsigned int body_id;
  hsize_t dimbody[1];
  H5LTget_attribute_uint(group_id,"body","ID",&body_id);

  H5LTget_dataset_info(group_id,"body", dimbody,NULL,NULL);
  double body_params[dimbody[0]];
  H5LTread_dataset_double(group_id,"body", body_params);

  hsize_t dimtrack[1];
  H5LTget_dataset_info(group_id,"track", dimtrack ,NULL,NULL);
  std::unique_ptr<double[]> track_params(new double[dimtrack[0]]);
  H5LTread_dataset_double(group_id,"track", track_params.get());

  double x_current;
  H5LTget_attribute_double(group_id,"track","X",&x_current);

  // setting body and track
  SetBodyTrack(body_id,dimbody[0],body_params,dimtrack[0],track_params.get());

  // set trayectory to current time
  track->SetX(x_current);

  // initializing nuSQUIDS
  if (ne == 1){
    if(not inusquids)
      init(squids_time_initial);
    Set_E(energy_data[0]);
  }
  else {
    init(E_range,false,squids_time_initial);
  }
  // reset current squids time
  Set_t(squids_time);
  // set time offset
  time_offset = squids_time - track->GetX();

  // evolve projectors to current time
  EvolveProjectors(squids_time);
  // reading state
  H5LTget_dataset_info(group_id,"neustate", dims,NULL,NULL);
  std::unique_ptr<double[]> neudata(new double[dims[0]*dims[1]]);
  H5LTread_dataset_double(group_id,"neustate", neudata.get());

  H5LTget_dataset_info(group_id,"aneustate", dims,NULL,NULL);
  std::unique_ptr<double[]> aneudata(new double[dims[0]*dims[1]]);
  H5LTread_dataset_double(group_id,"aneustate", aneudata.get());

  for(unsigned int ie = 0; ie < dims[0]; ie++){
    for (unsigned int j = 0; j < dims[1]; j ++){
      if (NT == neutrino)
        state[ie].rho[0][j] = neudata[ie*dims[1]+j];
      else if ( NT == antineutrino)
        state[ie].rho[0][j] = aneudata[ie*dims[1]+j];
      else if ( NT == both ){
        state[ie].rho[0][j] = neudata[ie*dims[1]+j];
        state[ie].rho[1][j] = aneudata[ie*dims[1]+j];
      }
    }
  }

  if(iinteraction){
    // if intereactions will be used then reading cross section information
    hid_t xs_grp;
    if ( cross_section_grp_loc == "") {
      xs_grp = H5Gopen(group_id, "crosssections", H5P_DEFAULT);
    } else {
      xs_grp = H5Gopen(root_id, cross_section_grp_loc.c_str(), H5P_DEFAULT);
    }

    // initialize vectors
    int_struct = std::make_shared<InteractionStructure>();
    InitializeInteractionVectors();

    // sigma_CC and sigma_NC

    hsize_t XSdim[3];
    H5LTget_dataset_info(xs_grp,"sigmacc", XSdim,NULL,NULL);

    std::unique_ptr<double[]> xsCC(new double[XSdim[0]*XSdim[1]*XSdim[2]]);
    H5LTread_dataset_double(xs_grp,"sigmacc", xsCC.get());
    std::unique_ptr<double[]> xsNC(new double[XSdim[0]*XSdim[1]*XSdim[2]]);
    H5LTread_dataset_double(xs_grp,"sigmanc", xsNC.get());

    for ( unsigned int rho = 0; rho < nrhos; rho ++){
      for ( unsigned int flv = 0; flv < numneu; flv ++){
          for ( unsigned int ie = 0; ie < ne; ie ++){
            int_struct->sigma_CC[rho][flv][ie] = xsCC[rho*(numneu*ne) +  flv*ne + ie];
            int_struct->sigma_NC[rho][flv][ie] = xsNC[rho*(numneu*ne) +  flv*ne + ie];
          }
      }
    }

    // dNdE_CC and dNdE_NC
    hsize_t dXSdim[4];
    H5LTget_dataset_info(xs_grp,"dNdEcc", dXSdim,NULL,NULL);

    std::unique_ptr<double[]> dxsCC(new double[dXSdim[0]*dXSdim[1]*dXSdim[2]*dXSdim[3]]);
    H5LTread_dataset_double(xs_grp,"dNdEcc", dxsCC.get());
    std::unique_ptr<double[]> dxsNC(new double[dXSdim[0]*dXSdim[1]*dXSdim[2]*dXSdim[3]]);
    H5LTread_dataset_double(xs_grp,"dNdEnc", dxsNC.get());

    for( unsigned int rho = 0; rho < nrhos; rho++){
      for( unsigned int flv = 0; flv < numneu; flv++){
          for( unsigned int e1 = 0; e1 < ne; e1++){
              for( unsigned int e2 = 0; e2 < e1; e2++){
                int_struct->dNdE_CC[rho][flv][e1][e2] = dxsCC[rho*(numneu*ne*ne) +  flv*ne*ne + e1*ne + e2];
                int_struct->dNdE_NC[rho][flv][e1][e2] = dxsNC[rho*(numneu*ne*ne) +  flv*ne*ne + e1*ne + e2];
              }
          }
      }
    }

    // sigma_GR
    {
      std::vector<double> xsGR(XSdim[2]);
      H5LTread_dataset_double(xs_grp,"sigmagr", xsGR.data());
      std::copy(xsGR.begin(), xsGR.end(), int_struct->sigma_GR.begin());
    }
    // dNdE_GR
    {
      std::vector<double> dxsGR(dXSdim[2]*dXSdim[3]);
      H5LTread_dataset_double(xs_grp,"dNdEgr", dxsGR.data());
      std::copy(dxsGR.begin(), dxsGR.end(), int_struct->dNdE_GR.begin());
    }

    // invlen_tau
    hsize_t iltdim[1];
    H5LTget_dataset_info(xs_grp,"invlentau", iltdim,NULL,NULL);
    std::unique_ptr<double[]> invlentau(new double[iltdim[0]]);
    H5LTread_dataset_double(xs_grp,"invlentau", invlentau.get());
    for(unsigned int ie = 0; ie < ne; ie ++)
      int_struct->invlen_tau[ie] = invlentau[ie];

    // dNdE_tau_all,dNdE_tau_lep
    hsize_t dNdEtaudim[2];
    H5LTget_dataset_info(xs_grp,"dNdEtauall", dNdEtaudim,NULL,NULL);
    
    std::unique_ptr<double[]> dNdEtauall(new double[dNdEtaudim[0]*dNdEtaudim[1]]);
    H5LTread_dataset_double(xs_grp,"dNdEtauall", dNdEtauall.get());
    std::unique_ptr<double[]> dNdEtaulep(new double[dNdEtaudim[0]*dNdEtaudim[1]]);
    H5LTread_dataset_double(xs_grp,"dNdEtaulep", dNdEtaulep.get());

    for( unsigned int e1 = 0; e1 < ne; e1++){
        for( unsigned int e2 = 0; e2 < e1; e2++){
          int_struct->dNdE_tau_all[e1][e2] = dNdEtauall[e1*ne + e2];
          int_struct->dNdE_tau_lep[e1][e2] = dNdEtaulep[e1*ne + e2];
        }
    }
  }

  // read from user parameters
  hid_t user_parameters_id = H5Gopen(group_id, "user_parameters", H5P_DEFAULT);
  //H5Eset_auto (H5E_DEFAULT,(H5E_auto_t) H5Eprint,stderr);
  AddToReadHDF5(user_parameters_id);
  //H5Eset_auto (H5E_DEFAULT,NULL, NULL);
  H5Gclose(user_parameters_id);

  // close HDF5 file
  H5Gclose ( group_id );
  // close root and file
  H5Gclose ( root_id );
  H5Fclose (file_id);

  // we assume that this was created with the writer and got to this point!
  istate = true;
  ienergy = true;
  itrack = true;
  ibody = true;
  // initialize H0
  iniH0();
}


void nuSQUIDS::SetBodyTrack(unsigned int body_id, unsigned int body_params_len, double body_params[], unsigned int track_params_len, double track_params[]){
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
          const unsigned int xn = body_params_len/3;
          std::vector<double> xx(xn),rho(xn),ye(xn);
          for(unsigned int i = 0; i < xn; i++){
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
          // track_param[2] corresponds to the zenith angle
          track = std::make_shared<EarthAtm::Track>(track_params[2]);
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

unsigned int nuSQUIDS::GetNumNeu() const{
  return numneu;
}

unsigned int nuSQUIDS::GetNumRho() const{
  return nrhos;
}

void nuSQUIDS::ProgressBar() const{
  double progress = (track->GetX()-track->GetInitialX())/(track->GetFinalX() - track->GetInitialX());
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
  if ( NT != both and opt )
    throw std::runtime_error("nuSQUIDS::Error::Cannot set TauRegeneration to True when NT != 'both'.");
  tauregeneration = opt;
}

void nuSQUIDS::Set_ProgressBar(bool opt){
  progressbar = opt;
}


void nuSQUIDS::Set_IncludeOscillations(bool opt){

  ioscillations = opt;
}

std::shared_ptr<Track> nuSQUIDS::GetTrack() const{
  return track;
}

std::shared_ptr<Body> nuSQUIDS::GetBody() const{
  return body;
}

void nuSQUIDS::Set_MixingAngle( unsigned int i, unsigned int j, double val){
  if ( i > numneu or j > numneu)
    throw std::invalid_argument("nuSQUIDS::Set_MixingAngle::Error: Mixing angle index greater than number of neutrino flavors.");
  params.SetMixingAngle(i,j,val);
  istate = false;
}

double nuSQUIDS::Get_MixingAngle( unsigned int i, unsigned int j) const {
  if ( i > numneu or j > numneu)
    throw std::invalid_argument("nuSQUIDS::Set_MixingAngle::Error: Mixing angle index greater than number of neutrino flavors.");
  return params.GetMixingAngle(i,j);
}

void nuSQUIDS::Set_CPPhase( unsigned int i, unsigned int j, double val){
  if ( i > numneu or j > numneu)
    throw std::invalid_argument("nuSQUIDS::Set_CPPhase::Error: CP phase index greater than number of neutrino flavors.");
  params.SetPhase(i,j,val);
  istate = false;
}

double nuSQUIDS::Get_CPPhase( unsigned int i, unsigned int j) const {
  if ( i > numneu or j > numneu)
    throw std::invalid_argument("nuSQUIDS::Set_CPPhase::Error: CP phase index greater than number of neutrino flavors.");
  return params.GetPhase(i,j);
}

void nuSQUIDS::Set_SquareMassDifference( unsigned int i, double val){
  if ( i > numneu )
    throw std::invalid_argument("nuSQUIDS::Set_SquareMassDifference::Error: Inder greater than number of neutrino flavors.");
  params.SetEnergyDifference(i,val);
  istate = false;
}

double nuSQUIDS::Get_SquareMassDifference( unsigned int i ) const {
  if ( i > numneu )
    throw std::invalid_argument("nuSQUIDS::Set_SquareMassDifference::Error: Inder greater than number of neutrino flavors.");
  return params.GetEnergyDifference(i);
}

void nuSQUIDS::Set_MixingParametersToDefault(void){
  // set parameters as in arXiv:1409.5439 NO
  // but with delta_CP = 0.0
  Set_MixingAngle(0,1,0.583996); // th12
  Set_MixingAngle(0,2,0.148190); // th13
  Set_MixingAngle(1,2,0.737324); // th23

  Set_SquareMassDifference(1,7.5e-05); // dm^2_21
  Set_SquareMassDifference(2,0.00257); // dm^2_31

  Set_CPPhase(0,2,0.0); // delta_13 = diract cp phase
}

void nuSQUIDS::Set_Basis(Basis b){
  if ( b == flavor )
    throw std::runtime_error("nuSQUIDS::Set_Basis::Error: solution basis can only be nuSQUIDS::mass or nuSQUIDS::interaction.");
  basis = b;
}

nuSQUIDS::~nuSQUIDS(){}

nuSQUIDS::nuSQUIDS(nuSQUIDS&& other):
squids::SQuIDS(std::move(other)),
basis(other.basis),
numneu(other.numneu),
ne(other.ne),
E_range(std::move(other.E_range)),
delE(std::move(other.delE)),
ncs(std::move(other.ncs)),
tdc(std::move(other.tdc)),
int_struct(std::move(other.int_struct)),
taubr_lep(other.taubr_lep),
tau_lifetime(other.tau_lifetime),
tau_mass(other.tau_mass),
tau_reg_scale(other.tau_reg_scale),
positivization_scale(other.positivization_scale),
body(other.body),
track(other.track),
DM2(other.DM2),
H0_array(std::move(other.H0_array)),
b0_proj(std::move(other.b0_proj)),
b1_proj(std::move(other.b1_proj)),
evol_b0_proj(std::move(other.evol_b0_proj)),
evol_b1_proj(std::move(other.evol_b1_proj)),
inusquids(other.inusquids),
ibody(other.ibody),
ienergy(other.ienergy),
itrack(other.itrack),
istate(other.istate),
iinteraction(other.iinteraction),
elogscale(other.elogscale),
tauregeneration(other.tauregeneration),
positivization(other.positivization),
progressbar(other.progressbar),
progressbar_count(other.progressbar_count),
progressbar_loop(other.progressbar_loop),
NT(other.NT)
{
  other.inusquids=false; //other is no longer usable, since we stole its contents
}

nuSQUIDS& nuSQUIDS::operator=(nuSQUIDS&& other){
  if(&other==this)
    return(*this);

  squids::SQuIDS::operator=(std::move(other));

  basis = other.basis;
  numneu = other.numneu;
  ne = other.ne;
  E_range = std::move(other.E_range);
  delE = std::move(other.delE);
  ncs = other.ncs;
  tdc = other.tdc;
  int_struct = std::move(other.int_struct);
  taubr_lep = other.taubr_lep;
  tau_lifetime = other.tau_lifetime;
  tau_mass = other.tau_mass;
  tau_reg_scale = other.tau_reg_scale;
  positivization_scale = other.positivization_scale;
  body = other.body;
  track = other.track;
  DM2 = other.DM2;
  H0_array = std::move(other.H0_array);
  b0_proj = std::move(other.b0_proj);
  b1_proj = std::move(other.b1_proj);
  evol_b0_proj = std::move(other.evol_b0_proj);
  evol_b1_proj = std::move(other.evol_b1_proj);

  inusquids = other.inusquids;
  ibody = other.ibody;
  ienergy = other.ienergy;
  itrack = other.itrack;
  istate = other.istate;
  iinteraction = other.iinteraction;
  elogscale = other.elogscale;
  tauregeneration = other.tauregeneration;
  positivization = other.positivization;
  progressbar = other.progressbar;
  progressbar_count = other.progressbar_count;
  progressbar_loop = other.progressbar_loop;

  NT = other.NT;

  // initial nusquids object render useless
  other.inusquids = false;

  return(*this);
}

} // close namespace
