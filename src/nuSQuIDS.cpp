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

#include <nuSQuIDS/nuSQuIDS.h>
#include <cstring>
#include <map>
#include <sstream>
#include <H5Dpublic.h>
#include <H5Epublic.h>
#include <H5Spublic.h>
#include "nuSQuIDS/resources.h"
#include "nuSQuIDS/version.h"

//Enable/disable finer-grained labels which can be useful for profiling
#define USE_DEBUG_LABELS 1
#if USE_DEBUG_LABELS
  #define NUSQUIDS_DEBUG_LABEL(tag) \
  asm("nusquids_debug_" #tag ":" :::)
#else
  #define NUSQUIDS_DEBUG_LABEL(tag)
#endif

namespace nusquids{


namespace {

std::string getStringAttribute(hid_t loc_id, const char* obj_name, const char* attr_name){
  herr_t err;
  int rank;
  err=H5LTget_attribute_ndims(loc_id,obj_name,attr_name,&rank);
  if(err<0)
    throw std::runtime_error("Failed to get dimensionality of attribute " + std::string(attr_name)+" of object "+std::string(obj_name));
  if(rank!=1)
    throw std::runtime_error("Attribute "+std::string(attr_name) +" of object "+std::string(obj_name)+" is not one-dimensional");
  hsize_t dims;
  H5T_class_t type_class;
  size_t type_size;
  err=H5LTget_attribute_info(loc_id,obj_name,attr_name,&dims,&type_class,&type_size);
  if(err<0)
    throw std::runtime_error("Failed to get info on attribute "+std::string(attr_name)+" of object "+std::string(obj_name));
  std::unique_ptr<char[]> buffer(new char[type_size+1]);
  if(!buffer)
    throw std::runtime_error("Failed to allocate space to read string atribute");
  err=H5LTget_attribute_string(loc_id,obj_name,attr_name,buffer.get());
  if(err<0)
    throw std::runtime_error("Failed to read attribute "+std::string(attr_name)+" of object "+std::string(obj_name));
  return(std::string(buffer.get(),type_size));
}

} // close unnamed namespace

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
  if ( numneu < 3)
    throw std::runtime_error("nuSQUIDS::Error::Minimum number of neutrinos is three");

  nsun = numneu;

  //initialize SQUIDS
  ini(ne,numneu,1,0,xini);
  Set_CoherentRhoTerms(true);
  Set_h(1.*params.km);
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

  //precompute this product for HI to avoid repeating expensive pow() calls.
  HI_constants = params.sqrt2*params.GF*params.Na*pow(params.cm,-3);

  //===============================
  // END                         //
  //===============================
  inusquids = true;
}

void nuSQUIDS::Set_E(double Enu){
  if( ne != 1 )
    throw std::runtime_error("nuSQUIDS::Error:Can only use Set_E in single energy mode.");
  E_range = marray<double,1>{1};
  E_range[0] = Enu;
  Set_xrange(std::vector<double>{Enu});

  // energy is initialized
  ienergy = true;
  // state is invalided, because hamiltonian changes.
  istate = false;
}

void nuSQUIDS::init(marray<double,1> E_vector, double xini){
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
  if ( numneu < 3)
    throw std::runtime_error("nuSQUIDS::Error::Minimum number of neutrinos is three");

  nsun = numneu;

  ne = E_vector.size();
  assert(ne>0);

  //===============================
  // BEGIN                       //
  //===============================

  try {
    // initialize SQUIDS
    ini(ne,numneu,nrhos,0,xini);
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    throw std::runtime_error("nuSQUIDS::init : Failed while trying to initialize SQuIDS.");
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

  positivization_scale = 300.0*params.km;

  if(iinteraction){
    Set_NonCoherentRhoTerms(true);
    Set_OtherRhoTerms(true);
  }
  
  //precompute this product for HI to avoid repeating expensive pow() calls.
  HI_constants = params.sqrt2*params.GF*params.Na*pow(params.cm,-3);

  //===============================
  // END                         //
  //===============================
}

void nuSQUIDS::InitializeInteractionVectors(unsigned int nTargets){
  //Note that by rounding up the size of each array's final dimension we ensure
  //that if the beginning of the array is aligned, each segment starting with an
  //entry whose final index is zero will be as well.
  
  // initialize cross section and interaction arrays
  size_t rounded_ne=round_up_to_aligned(ne);
  int_struct->dNdE_NC.resize(std::vector<size_t>{nTargets,nrhos,numneu,ne,rounded_ne});
  int_struct->dNdE_CC.resize(std::vector<size_t>{nTargets,nrhos,numneu,ne,rounded_ne});
  int_struct->dNdE_GR.resize(std::vector<size_t>{ne,rounded_ne});
  // initialize cross section arrays
  int_struct->sigma_CC.resize(std::vector<size_t>{nTargets,nrhos,numneu,rounded_ne});
  int_struct->sigma_NC.resize(std::vector<size_t>{nTargets,nrhos,numneu,rounded_ne});
  int_struct->sigma_GR.resize(std::vector<size_t>{rounded_ne});
  // initialize the tau decay and interaction array
  int_struct->dNdE_tau_all.resize(std::vector<size_t>{nrhos,ne,rounded_ne});
  int_struct->dNdE_tau_lep.resize(std::vector<size_t>{nrhos,ne,rounded_ne});
}

void nuSQUIDS::PreDerive(double x){
  track->SetX(x-time_offset);
  current_ye = body->ye(*track);
  current_density = body->density(*track);
  if(enable_neutrino_sources){
    body->injected_neutrino_flux(current_external_flux,*track,*this);
    assert(current_external_flux.extent(0) == ne && current_external_flux.extent(1) == nrhos && current_external_flux.extent(2) == numneu);
  }
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
  std::unique_ptr<double[]> evol_buf(new double[H0_array[0].GetEvolveBufferSize()]);
  for(unsigned int ei = 0; ei < ne; ei++){
    H0_array[ei].PrepareEvolve(evol_buf.get(),x-Get_t_initial());
    if (evol_lowpass_cutoff > 0){
      // std::cout << "evol_buf: ";
      // for (unsigned int i = 0; i < H0_array[ei].GetEvolveBufferSize(); i++){
      //   std::cout << evol_buf.get()[i] << "  ";
      // }
      // std::cout << std::endl;
      H0_array[ei].LowPassFilter(evol_buf.get(), evol_lowpass_cutoff, evol_lowpass_scale);
      // std::cout << "evol_buf after filter: ";
      // for (unsigned int i = 0; i < H0_array[ei].GetEvolveBufferSize(); i++){
      //   std::cout << evol_buf.get()[i] << "  ";
      // }
      // std::cout << std::endl;
    }
    for(unsigned int rho = 0; rho < nrhos; rho++){
      for(unsigned int flv = 0; flv < numneu; flv++){
        // will only evolve the flavor projectors
        evol_b1_proj[rho][flv][ei] = squids::detail::guarantee
          <squids::detail::NoAlias | squids::detail::EqualSizes>
          (b1_proj[rho][flv].Evolve(evol_buf.get()));
      }
    }
  }
  return;
}

squids::SU_vector nuSQUIDS::H0(double Enu, unsigned int irho) const{
  if (not ioscillations)
    return squids::SU_vector(nsun);
  return DM2*(0.5/Enu);
}

squids::SU_vector nuSQUIDS::HI(unsigned int ie, unsigned int irho) const{
    double CC = HI_constants*current_density*current_ye;
    double NC;

    if (current_ye < 1.0e-10)
      NC = HI_constants*current_density;
    else
      NC = CC*(-0.5*(1.0-current_ye)/current_ye);
    // antineutrino matter potential flips sign
    if ((irho == 1 and NT==both) or NT==antineutrino){
      CC*=-1;
      NC*=-1;
    }

    // construct potential in flavor basis
    squids::SU_vector potential((CC+NC)*evol_b1_proj[irho][0][ie]);
    potential += squids::detail::guarantee
      <squids::detail::NoAlias | squids::detail::EqualSizes | squids::detail::AlignedStorage>
      (NC*evol_b1_proj[irho][1][ie]);
    potential += squids::detail::guarantee
      <squids::detail::NoAlias | squids::detail::EqualSizes | squids::detail::AlignedStorage>
      (NC*evol_b1_proj[irho][2][ie]);

    if (basis == mass)
      potential += H0_array[ie];

    return potential;
}

squids::SU_vector nuSQUIDS::GammaRho(unsigned int ei,unsigned int index_rho) const{
    if (not iinteraction)
      return squids::SU_vector(nsun);

    squids::SU_vector V(evol_b1_proj[index_rho][0][ei]*(0.5*int_state.invlen_INT[index_rho][0][ei]));
    V += evol_b1_proj[index_rho][1][ei]*(0.5*int_state.invlen_INT[index_rho][1][ei]);
    V += evol_b1_proj[index_rho][2][ei]*(0.5*int_state.invlen_INT[index_rho][2][ei]);

    return V;
}

squids::SU_vector nuSQUIDS::InteractionsRho(unsigned int e1,unsigned int index_rho) const{
  if(iinteraction && !ioscillations) //can use precomputed result from UpdateInteractions
    return(interaction_cache[index_rho][e1]);

  if (not iinteraction and not enable_neutrino_sources)
    return squids::SU_vector(nsun);
  else if (not iinteraction) {
    squids::SU_vector external_flux_contribution(current_external_flux[e1][index_rho][0]*evol_b1_proj[index_rho][0][e1]);
    for(unsigned int i = 1; i < numneu; i++)
      external_flux_contribution+=squids::detail::guarantee
                                  <squids::detail::NoAlias | squids::detail::EqualSizes | squids::detail::AlignedStorage>
                                  (current_external_flux[e1][index_rho][i]*evol_b1_proj[index_rho][i][e1]);
    return external_flux_contribution;
  }

  // NC interactions
  double factor_e  =nc_factors[index_rho][0][e1];
  double factor_mu =nc_factors[index_rho][1][e1];
  double factor_tau=nc_factors[index_rho][2][e1];
  // Tau regeneration
  if(tauregeneration){
    factor_e  +=   tau_lep_decays[index_rho][e1];
    factor_mu +=   tau_lep_decays[index_rho][e1];
    factor_tau+=tau_hadlep_decays[index_rho][e1];
  }
  // Glashow resonance for electron antineutrinos
  if (iglashow && ((NT == both and index_rho == 1) or NT == antineutrino)) {
    factor_e  +=gr_factors[e1];
    factor_mu +=gr_factors[e1];
    factor_tau+=gr_factors[e1];
  }

  // add explicitly
  if(enable_neutrino_sources){
    factor_e  +=current_external_flux[e1][index_rho][0];
    factor_mu +=current_external_flux[e1][index_rho][1];
    factor_tau+=current_external_flux[e1][index_rho][2];
  }

  // Add the weighted projectors
  squids::SU_vector interaction_term(factor_e*evol_b1_proj[index_rho][0][e1]);
  interaction_term+=squids::detail::guarantee
                    <squids::detail::NoAlias | squids::detail::EqualSizes | squids::detail::AlignedStorage>
                    (factor_mu*evol_b1_proj[index_rho][1][e1]);
  interaction_term+=squids::detail::guarantee
                    <squids::detail::NoAlias | squids::detail::EqualSizes | squids::detail::AlignedStorage>
                    (factor_tau*evol_b1_proj[index_rho][2][e1]);

  if(enable_neutrino_sources){
    for(unsigned int i = 3; i < numneu; i++){
      interaction_term+=squids::detail::guarantee
                        <squids::detail::NoAlias | squids::detail::EqualSizes | squids::detail::AlignedStorage>
                        (current_external_flux[e1][index_rho][i]*evol_b1_proj[index_rho][i][e1]);
    }
  }

  return interaction_term;
}

double nuSQUIDS::GammaScalar(unsigned int ei, unsigned int iscalar) const{
  return 0.0;
}

double nuSQUIDS::InteractionsScalar(unsigned int ei, unsigned int iscalar) const{
  return 0.;
}
  
std::vector<double> nuSQUIDS::GetTargetNumberFractions() const{
  if(int_struct->targets.size()==1)
    return {1}; //with only one target type, it must make up all of the material
  //note that here we assume the order of targets defined in InitializeInteractions is proton, neutron
  return {current_ye,1-current_ye};
}
  
std::vector<double> nuSQUIDS::GetTargetNumberDensities() const{
  //start from the mass density, converted to natural units
  double density = (params.gr*pow(params.cm,-3))*current_density;
  
  //check whether we are treating the medium as isoscalar
  if(int_struct->targets.size()==1 && int_struct->targets[0]==isoscalar_nucleon){
    double num_nuc = density*2.0/(params.electron_mass+params.proton_mass+params.neutron_mass);
    if(num_nuc < 1.0e-10 )
      num_nuc = params.Na*pow(params.cm,-3)*1.0e-10;
    return {num_nuc};
  }
  
  //otherwise, treat protons and neutrons separately
  //the electron fraction is electrons/nucleon:
  //ye = ne / (np + nn)
  //The mass density is:
  //density = Me * ne + Mp * np + Mn * nn
  //assuming neutral matter ne == np
  //ye = np / (np + nn)
  //density = (Me + Mp) * np + Mn * nn
  //
  //(np + nn) * ye = np
  //ye * np + ye * nn = np
  //ye * nn = (1 - ye) * np
  //in the case that ye != 0:
  //nn = (1 - ye) * np / ye
  //
  //density = (Me + Mp) * np + Mn * (1 - ye) * np / ye
  //density = ((Me + Mp) + Mn * (1 - ye) / ye) * np
  //np = density / ((Me + Mp) + Mn * (1 - ye) / ye)
  
  double ye = current_ye;
  if(ye == 0){
    //apparently the medium is all neutrons
    double nn=density/params.neutron_mass;
    return {0,nn};
  }
  double np = density / (params.electron_mass + params.proton_mass + params.neutron_mass*((1-ye)/ye));
  double nn = (1-ye)*np/ye;
  //note that here we assume the order of targets defined in InitializeInteractions
  return {np,nn};
}
  
void nuSQUIDS::SetUpInteractionCache(){
  //check whether initialization has already been done
  if(interaction_cache.extent(0)==nrhos &&
     interaction_cache.extent(1)==ne &&
     interaction_cache_store)
    return;
  const size_t nvectors=nrhos*ne;
  const size_t vector_size=numneu*numneu;
  //std::cout << "N Vectors: " << nvectors << std::endl;
  //std::cout << "Vector size: " << vector_size << std::endl;
  
  //Both of these variables have units of sizeof(double), not bytes!
  size_t headroom=0;
  size_t padding=0;
  
  //If doubles are already at least 32 byte aligned, we don't need to do
  //anything, and in the SU(1) case we don't care about alignment.
  if(alignof(double)<32 && numneu>1){
    assert(sizeof(double)<32); //we're going to assume this throughout
  //We want a very particular alignment: if our vectors have more than one
  //component (SU(2) or greater) we want (1) the components to be 32 byte
  //aligned if the number of components is even, or (2) if the number of
  //components is odd we want the _second_ component to be 32 byte aligned.
  //We assume that operator new[] will return a block which is aligned to
  //sizeof(double), but not necessarily moreso.
  //In both cases we may need to skip up to 32-sizeof(double) bytes at the
  //start of the block.
  //Next, we may need padding between the vectors' individual chunks. In case
  //(1) this will be just
  //((vector_size*sizeof(double))%32 ? 32-((vector_size*sizeof(double))%32) : 0)
  //bytes and in case
  //(2) this will be 32-sizeof(double)-((vector_size-1)*sizeof(double))%32,
  //plus 32 if negative.
  
  headroom = (32-sizeof(double))/sizeof(double);
  if(vector_size%2){ //odd number of elements; case (2)
    size_t overhang=((vector_size-1)*sizeof(double))%32;
    size_t gap=32-overhang;
    //std::cout << " overhang: " << overhang << std::endl;
    //std::cout << " gap: " << gap << std::endl;
    if(gap==sizeof(double)) //the first component of the vector will exactly fit
      padding=0;
    else if(gap<sizeof(double)){ //need to make room for first
      overhang=sizeof(double)-gap;
      padding=(32-overhang)/sizeof(double);
    }
    else if(gap>sizeof(double)){
      padding=(gap-sizeof(double))/sizeof(double);
    }
  }
  else{ //even number; case (1)
    size_t overhang=(vector_size*sizeof(double))%32;
    if(overhang)
      padding=(32-overhang)/sizeof(double);
  }
  }
  
  //std::cout << "headroom needed: " << headroom << std::endl;
  //std::cout << "padding needed: " << padding << std::endl;
  
  //headroom at the front, and room for nvectors vectors, with (nvectors-1)
  //sections of padding between them.
  interaction_cache_store_size=headroom + nvectors*vector_size + (nvectors-1)*padding;
  //std::cout << "Allocating " << interaction_cache_store_size*sizeof(double) << " bytes" << std::endl;
  //std::cout << " space wasted achieving alignment: " << (headroom + (nvectors-1)*padding)*sizeof(double) << " bytes" << std::endl;
  interaction_cache_store.reset(new double[interaction_cache_store_size]);
  
  interaction_cache.resize(std::initializer_list<size_t>{nrhos,ne});
  //std::cout << "Initial pointer is " << interaction_cache_store.get() << std::endl;
  double* ptr=interaction_cache_store.get();
  if(alignof(double)<32 && numneu>1){
    if(vector_size%2){ //odd number of elements; case (2)
      if((intptr_t)(ptr+1)%32)
        ptr+=(32-((intptr_t)(ptr+1)%32))/sizeof(double);
    }
    else if((intptr_t)ptr%32) //even number of elements; case (1)
      ptr+=(32-((intptr_t)ptr%32))/sizeof(double);
  }
  //std::cout << "Pointer after adjustment is " << ptr << std::endl;
  for(size_t i=0; i<nrhos; i++){
    for(size_t j=0; j<ne; j++){
      if(alignof(double)<32 && numneu>1){
        size_t alignment=(intptr_t)ptr%32;
        if((numneu%2==1 && alignment!=32-sizeof(double)) || (numneu%2==0 && alignment!=0)){
          //std::cout << "Pointer is " << ptr << ", alignment is " << alignment << std::endl;
          throw std::logic_error("Fatal error: Logic in nuSQUIDS::SetUpInteractionCache has failed to ensure 32 byte aligned pointers");
        }
      }
      interaction_cache[i][j]=squids::SU_vector(numneu,ptr);
      ptr+=vector_size+padding;
    }
  }
}

void nuSQUIDS::UpdateInteractions(){
  //make sure vectors are the right size
  size_t rounded_ne=round_up_to_aligned(ne);
  if(int_state.invlen_NC.extent(0)!=nrhos ||
     int_state.invlen_NC.extent(1)!=numneu ||
     int_state.invlen_NC.extent(2)!=rounded_ne){
    int_state.invlen_NC.resize(std::vector<size_t>{nrhos,numneu,rounded_ne});
    int_state.invlen_CC.resize(std::vector<size_t>{nrhos,numneu,rounded_ne});
    int_state.invlen_GR.resize(std::vector<size_t>{rounded_ne});
    int_state.invlen_INT.resize(std::vector<size_t>{nrhos,numneu,rounded_ne});
  }
  
  std::vector<double> targetDensities = GetTargetNumberDensities();
  std::vector<double> targetFractions = GetTargetNumberFractions();
  const unsigned int nTargetTypes = targetDensities.size();
  
  //A helper class for accumulating a sum into a buffer without needing a first pass 
  //to fill the buffer with zeros
  struct LinAcc{
    bool add;
    LinAcc(unsigned int iteration):add(iteration){}
    void operator()(double& target, double contribution) const{
      if(add)
        target+=contribution;
      else
        target=contribution;
    }
  };
  
  NUSQUIDS_DEBUG_LABEL(setting_invlens);
  for(unsigned int rho = 0; rho < nrhos; rho++){
    for(unsigned int flv = 0; flv < numneu; flv++){
      for(unsigned int e1 = 0; e1 < ne; e1++){
        for(unsigned int trg=0; trg < nTargetTypes; trg++){
          LinAcc acc(trg);
          acc(int_state.invlen_NC[rho][flv][e1], int_struct->sigma_NC[trg][rho][flv][e1]*targetDensities[trg]);
          acc(int_state.invlen_CC[rho][flv][e1], int_struct->sigma_CC[trg][rho][flv][e1]*targetDensities[trg]);
        }
        int_state.invlen_INT[rho][flv][e1] = int_state.invlen_NC[rho][flv][e1] + int_state.invlen_CC[rho][flv][e1];
      }
    }
  }
    
  // Add Glashow resonance if antineutrinos are in the mix
  if (iglashow && (NT == both or NT == antineutrino)) {
    NUSQUIDS_DEBUG_LABEL(setting_invlen_GR);
    assert(numneu > 0);
    unsigned int rho = (NT == both) ? 1 : 0;
    //if we are dealing in generic 'nucleons' due to isoscalar cross sections, scale by electron 
    //fraction to get electron number density. Otherwise, just take the proton number density. 
    double num_e = nTargetTypes==1?targetDensities[0]*current_ye:targetDensities[0];
    double* sigma_GR_ptr=&int_struct->sigma_GR[0];
    double* invlen_GR_ptr=&int_state.invlen_GR[0];
    double* invlen_INT_ptr=&int_state.invlen_INT[rho][0][0];
    SQUIDS_POINTER_IS_ALIGNED(sigma_GR_ptr,preferred_alignment*sizeof(double));
    SQUIDS_POINTER_IS_ALIGNED(invlen_GR_ptr,preferred_alignment*sizeof(double));
    SQUIDS_POINTER_IS_ALIGNED(invlen_INT_ptr,preferred_alignment*sizeof(double));
    for(unsigned int e1 = 0; e1 < ne; e1++){
      invlen_GR_ptr[e1] = sigma_GR_ptr[e1]*num_e;
      invlen_INT_ptr[e1] += invlen_GR_ptr[e1];
    }
  }

  //Create a local buffer, in the form of an std::vector holding values of type
  //buf_type, with space for at least buf_size elements. Align data storage
  //according to preferred_alignment, rounding up the actual size of the buffer
  //if necessary. All data is initialized to zero.
  //Produces a pointer named buf_name to the resulting underlying storage, and
  //hints to the compiler that this pointer points to aligned data.
  #define ALIGNED_LOCAL_BUFFER(buf_name,buf_type,buf_size) \
    std::vector<buf_type,aligned_allocator<buf_type>> buf_name ## _vec( \
      round_up_to_aligned(buf_size),0, \
      aligned_allocator<buf_type>(log2(preferred_alignment*sizeof(buf_type)))); \
    buf_type* buf_name=(buf_name ## _vec).data(); \
    SQUIDS_POINTER_IS_ALIGNED(buf_name,preferred_alignment*sizeof(buf_type)) \
  // ^ Note lack of trailing semicolon
  
  //Without oscillations, the entries in evol_b1_proj do not depend on energy.
  //We can exploit this to precalculate the information needed by InteractionsRho
  //performing operations on SU_vectors a number of times proportional to the number
  //of energies, rather than proportinal to the square.
  if(!ioscillations){
    auto rounded_ne=round_up_to_aligned(ne);
    marray<double,3,aligned_allocator<double>> flavor_factors({nrhos,3,rounded_ne},aligned_allocator<double>(log2(preferred_alignment*sizeof(double))));
    
    // NC interactions
    NUSQUIDS_DEBUG_LABEL(NC_contribution_no_osc);
    for(unsigned int rho = 0; rho < nrhos; rho++){
      //for each flavor
      for(unsigned int alpha_active : {0,1,2}){
        //accumulate the contribution of each energy e2 to each lower energy
        squids::SU_vector& projector=evol_b1_proj[rho][alpha_active][0];
        
        for(unsigned int e2=1; e2<ne; e2++){
          //the flux of the current flavor at e2
          double flux_a_e2=projector*estate[e2].rho[rho];
          //premultiply factors which do not depend on the lower energy e1
          flux_a_e2*=int_state.invlen_NC[rho][alpha_active][e2]*delE[e2-1];
          
          for(unsigned int trg=0; trg < nTargetTypes; trg++){
            double tFrac=targetFractions[trg];
            //the target fraction factor is logically more a part of dN/dE, but since it is common 
            //to all terms it is more efficient if we attach it to the flux
            double flux_a_e2_t=tFrac*flux_a_e2;
            double* factors=&flavor_factors[rho][alpha_active][0];
            SQUIDS_POINTER_IS_ALIGNED(factors,preferred_alignment*sizeof(double));
            double* dNdE_ptr=&int_struct->dNdE_NC[trg][rho][alpha_active][e2][0];
            SQUIDS_POINTER_IS_ALIGNED(dNdE_ptr,preferred_alignment*sizeof(double));
            for(unsigned int e1=0; e1<e2; e1++, dNdE_ptr++)
              factors[e1]+=flux_a_e2_t*(*dNdE_ptr);
          }
        }
      }
    }

    NUSQUIDS_DEBUG_LABEL(TR_contribution_no_osc);
    if(tauregeneration){
      assert(numneu >= 3);
      const unsigned int tau_flavor = 2;
      
      //first accumulate the flux of taus produced by interaction
      ALIGNED_LOCAL_BUFFER(    tau_decay_fluxes,double,ne);
      ALIGNED_LOCAL_BUFFER(tau_bar_decay_fluxes,double,ne);
      
      //without osciallations we can always just use the first projectors
      squids::SU_vector& projector_tau     = evol_b1_proj[0][tau_flavor][0];
      squids::SU_vector& projector_tau_bar = evol_b1_proj[1][tau_flavor][0];
      
      for(unsigned int en=1; en<ne; en++){ // loop over initial tau neutrino energies
        double nu_tau_flux     =     projector_tau*estate[en].rho[0];
        double nu_tau_bar_flux = projector_tau_bar*estate[en].rho[1];
        if(nu_tau_flux<=0 && nu_tau_bar_flux<=0)
          continue;
        double dEn = delE[en-1];
        double invlen_CC_tau     = int_state.invlen_CC[0][tau_flavor][en];
        double invlen_CC_tau_bar = int_state.invlen_CC[1][tau_flavor][en];
        for(unsigned int trg=0; trg < nTargetTypes; trg++){ //loop over target types
          //note that the flux is scaled by fraction of total targets
          double nu_tau_flux_invlen_CC     =     nu_tau_flux*    invlen_CC_tau*dEn*targetFractions[trg];
          double nu_tau_bar_flux_invlen_CC = nu_tau_bar_flux*invlen_CC_tau_bar*dEn*targetFractions[trg];
          double* dNdE_CC_tau     = &int_struct->dNdE_CC[trg][0][tau_flavor][en][0];
          double* dNdE_CC_tau_bar = &int_struct->dNdE_CC[trg][1][tau_flavor][en][0];
          SQUIDS_POINTER_IS_ALIGNED(dNdE_CC_tau    ,preferred_alignment*sizeof(double));
          SQUIDS_POINTER_IS_ALIGNED(dNdE_CC_tau_bar,preferred_alignment*sizeof(double));
          for(unsigned int et=1; et<en; et++){ // loop over intermediate tau energies
            double dEt = delE[et-1];
            tau_decay_fluxes[et]    +=    nu_tau_flux_invlen_CC*    dNdE_CC_tau[et]*dEt;
            tau_bar_decay_fluxes[et]+=nu_tau_bar_flux_invlen_CC*dNdE_CC_tau_bar[et]*dEt;
          }
        }
      }
      
      //then accumulate the contributions back to the neutrino fluxes from the taus decaying
      ALIGNED_LOCAL_BUFFER(    tau_lep_decays,double,ne);
      ALIGNED_LOCAL_BUFFER(tau_bar_lep_decays,double,ne);
      double* factors_nu_e      =&flavor_factors[0][0][0];
      double* factors_nu_e_bar  =&flavor_factors[1][0][0];
      double* factors_nu_mu     =&flavor_factors[0][1][0];
      double* factors_nu_mu_bar =&flavor_factors[1][1][0];
      double* factors_nu_tau    =&flavor_factors[0][2][0];
      double* factors_nu_tau_bar=&flavor_factors[1][2][0];
      SQUIDS_POINTER_IS_ALIGNED(factors_nu_tau,preferred_alignment*sizeof(double));
      SQUIDS_POINTER_IS_ALIGNED(factors_nu_tau_bar,preferred_alignment*sizeof(double));
      for(unsigned int et=1; et<ne; et++){ // loop over intermediate tau energies
        double* tau_all_ptr=&int_struct->dNdE_tau_all[0][et][0];
        double* tau_lep_ptr=&int_struct->dNdE_tau_lep[0][et][0];
        double* taubar_all_ptr=&int_struct->dNdE_tau_all[1][et][0];
        double* taubar_lep_ptr=&int_struct->dNdE_tau_lep[1][et][0];
        SQUIDS_POINTER_IS_ALIGNED(tau_all_ptr,preferred_alignment*sizeof(double));
        SQUIDS_POINTER_IS_ALIGNED(tau_lep_ptr,preferred_alignment*sizeof(double));
        SQUIDS_POINTER_IS_ALIGNED(taubar_all_ptr,preferred_alignment*sizeof(double));
        SQUIDS_POINTER_IS_ALIGNED(taubar_lep_ptr,preferred_alignment*sizeof(double));
#ifdef __clang__ //clang needs a little help with this one
//this pragma is only available in new enough clang versions, unfortunately, 
//Apple screws up the version numbers in their copies, so this is rather messy
#if (defined(__apple_build_version__) && __clang_major__>=7) || \
    (!defined(__apple_build_version__) && (__clang_major__>3 || (__clang_major__==3 && __clang_minor__>=7)))
        #pragma clang loop vectorize(assume_safety)
#endif
#endif
        for(unsigned int e1=0; e1<et; e1++){ // loop over final neutrino energies
          //Defer adding the leptonic decays to both the electron and muon sums
          //to reduce the amount of data written in this inner loop from ~2*ne^2
          //to ne^2 + 2*ne. Since the hadlep contributions only go into the tau
          //components we might as well put them there directly.
          tau_lep_decays[e1]     += tau_bar_decay_fluxes[et]*taubar_lep_ptr[e1];
          tau_bar_lep_decays[e1] +=     tau_decay_fluxes[et]*tau_lep_ptr[e1];
          factors_nu_tau[e1]     +=     tau_decay_fluxes[et]*tau_all_ptr[e1];
          factors_nu_tau_bar[e1] += tau_bar_decay_fluxes[et]*taubar_all_ptr[e1];
        }
      }
      for(unsigned int e1=0; e1<rounded_ne; e1++){ // loop over final neutrino energies
        double tau_lep_decays_e1       =       tau_lep_decays[e1];
        double tau_bar_lep_decays_e1   =   tau_bar_lep_decays[e1];
        factors_nu_e      [e1]+=       tau_lep_decays_e1;
        factors_nu_e_bar  [e1]+=   tau_bar_lep_decays_e1;
        factors_nu_mu     [e1]+=       tau_lep_decays_e1;
        factors_nu_mu_bar [e1]+=   tau_bar_lep_decays_e1;
      }
    } // end of TR handling
    
    NUSQUIDS_DEBUG_LABEL(GR_contribution_no_osc);
    if(iglashow && (NT == both or NT == antineutrino)){
      unsigned int rho = (NT==both) ? 1 : 0;
      squids::SU_vector projector_sum(nsun);
      squids::SU_vector projector_e = evol_b1_proj[rho][0][0];
      projector_sum += projector_e;
      projector_sum += evol_b1_proj[rho][1][0];
      projector_sum += evol_b1_proj[rho][2][0];
      ALIGNED_LOCAL_BUFFER(gr_factors,double,ne);
      for(unsigned int e2=1; e2<ne; e2++){
        double flux=projector_e*estate[e2].rho[rho];
        double flux_invlen_en=flux*int_state.invlen_GR[e2]*delE[e2-1];
        double* dNdE_GR_ptr=&int_struct->dNdE_GR[e2][0];
        SQUIDS_POINTER_IS_ALIGNED(dNdE_GR_ptr,preferred_alignment*sizeof(double));
        for(unsigned int e1=0; e1<e2; e1++)
          gr_factors[e1] += flux_invlen_en*dNdE_GR_ptr[e1];
      }
      double* factors_nu_e_bar  =&flavor_factors[rho][0][0];
      double* factors_nu_mu_bar =&flavor_factors[rho][1][0];
      double* factors_nu_tau_bar=&flavor_factors[rho][2][0];
      SQUIDS_POINTER_IS_ALIGNED(factors_nu_e_bar  ,preferred_alignment*sizeof(double));
      SQUIDS_POINTER_IS_ALIGNED(factors_nu_mu_bar ,preferred_alignment*sizeof(double));
      SQUIDS_POINTER_IS_ALIGNED(factors_nu_tau_bar,preferred_alignment*sizeof(double));
      for(unsigned int e1=0; e1<ne; e1++){
        double gr_factor=gr_factors[e1];
        factors_nu_e_bar  [e1]+=gr_factor;
        factors_nu_mu_bar [e1]+=gr_factor;
        factors_nu_tau_bar[e1]+=gr_factor;
      }
    }
    
    // Reset to zero
    memset(interaction_cache_store.get(),0,interaction_cache_store_size*sizeof(double));
    // Add the weighted projectors
    for(unsigned int rho = 0; rho < nrhos; rho++){
      NUSQUIDS_DEBUG_LABEL(summing_projectors_no_osc);
      squids::SU_vector projector_e   = evol_b1_proj[rho][0][0];
      squids::SU_vector projector_mu  = evol_b1_proj[rho][1][0];
      squids::SU_vector projector_tau = evol_b1_proj[rho][2][0];
      double* factors_nu_e  =&flavor_factors[rho][0][0];
      double* factors_nu_mu =&flavor_factors[rho][1][0];
      double* factors_nu_tau=&flavor_factors[rho][2][0];
      SQUIDS_POINTER_IS_ALIGNED(factors_nu_e  ,preferred_alignment*sizeof(double));
      SQUIDS_POINTER_IS_ALIGNED(factors_nu_mu ,preferred_alignment*sizeof(double));
      SQUIDS_POINTER_IS_ALIGNED(factors_nu_tau,preferred_alignment*sizeof(double));
      for(unsigned int e1=0; e1<ne; e1++){
        interaction_cache[rho][e1]+=squids::detail::guarantee
                                    <squids::detail::EqualSizes | squids::detail::AlignedStorage>
                                    (factors_nu_e[e1]  *projector_e);
        interaction_cache[rho][e1]+=squids::detail::guarantee
                                    <squids::detail::EqualSizes | squids::detail::AlignedStorage>
                                    (factors_nu_mu[e1] *projector_mu);
        interaction_cache[rho][e1]+=squids::detail::guarantee
                                    <squids::detail::EqualSizes | squids::detail::AlignedStorage>
                                    (factors_nu_tau[e1]*projector_tau);
      }

    if(enable_neutrino_sources){
      for(unsigned int rho = 0; rho < nrhos; rho++){
        for(unsigned int e1=0; e1<ne; e1++){
          for(unsigned int i = 0; i < numneu; i++){
            interaction_cache[rho][e1]+=squids::detail::guarantee
                                        <squids::detail::NoAlias | squids::detail::EqualSizes | squids::detail::AlignedStorage>
                                        (current_external_flux[e1][rho][i]*evol_b1_proj[rho][i][e1]);
          }
        }
      }
    }

    }

  } //end of no oscillation handling
  else{ //oscillations
    NUSQUIDS_DEBUG_LABEL(NC_contribution_osc);
    
    std::fill(nc_factors.begin(),nc_factors.end(),0.);
    for(unsigned int rho = 0; rho < nrhos; rho++){
      //for each flavor
      for(unsigned int alpha_active : {0,1,2}){
        //accumulate the contribution of each energy e2 to each lower energy
        for(unsigned int e2=1; e2<ne; e2++){
          //the flux of the current flavor at e2
          double flux_a_e2=evol_b1_proj[rho][alpha_active][e2]*estate[e2].rho[rho];
          //premultiply factors which do not depend on the lower energy e1
          flux_a_e2*=int_state.invlen_NC[rho][alpha_active][e2]*delE[e2-1];
          //for each type of target within the material
          for(unsigned int trg=0; trg<nTargetTypes; trg++){
            //scale by the fraction of the material which is this target type
            double tFrac=targetFractions[trg];
            double flux_a_e2_t=tFrac*flux_a_e2;
            double* dNdE_ptr=&int_struct->dNdE_NC[trg][rho][alpha_active][e2][0];
            SQUIDS_POINTER_IS_ALIGNED(dNdE_ptr,preferred_alignment*sizeof(double));
            double* nc_factor_ptr=&nc_factors[rho][alpha_active][0];
            SQUIDS_POINTER_IS_ALIGNED(nc_factor_ptr,preferred_alignment*sizeof(double));
            for(unsigned int e1=0; e1<e2; e1++, dNdE_ptr++, nc_factor_ptr++)
              *nc_factor_ptr+=flux_a_e2_t*(*dNdE_ptr);
          } 
        }
      }
    }
    
    NUSQUIDS_DEBUG_LABEL(TR_contribution_osc);
    if(tauregeneration){
      assert(numneu >= 3);
      const unsigned int tau_flavor = 2;
      //first accumulate the flux of taus produced by interaction
      ALIGNED_LOCAL_BUFFER(tau_decay_fluxes,double,ne);
      ALIGNED_LOCAL_BUFFER(tau_bar_decay_fluxes,double,ne);
      //we must use the actual projector for each energy, no shortcuts!
      for(unsigned int en=1; en<ne; en++){ // loop over initial tau neutrino energies
        double nu_tau_flux     = evol_b1_proj[0][tau_flavor][en]*estate[en].rho[0];
        double nu_tau_bar_flux = evol_b1_proj[1][tau_flavor][en]*estate[en].rho[1];
        if(nu_tau_flux<=0 && nu_tau_bar_flux<=0)
          continue;
        double dEn = delE[en-1];
        double invlen_CC_tau     = int_state.invlen_CC[0][tau_flavor][en];
        double invlen_CC_tau_bar = int_state.invlen_CC[1][tau_flavor][en];
        for(unsigned int trg=0; trg < nTargetTypes; trg++){
          //for convenience scale fluxes by target fraction, as usual
          double nu_tau_flux_invlen_CC     =     nu_tau_flux*    invlen_CC_tau*dEn*targetFractions[trg];
          double nu_tau_bar_flux_invlen_CC = nu_tau_bar_flux*invlen_CC_tau_bar*dEn*targetFractions[trg];
          double* dNdE_CC_tau     = &int_struct->dNdE_CC[trg][0][tau_flavor][en][0];
          double* dNdE_CC_tau_bar = &int_struct->dNdE_CC[trg][1][tau_flavor][en][0];
          SQUIDS_POINTER_IS_ALIGNED(dNdE_CC_tau    ,preferred_alignment*sizeof(double));
          SQUIDS_POINTER_IS_ALIGNED(dNdE_CC_tau_bar,preferred_alignment*sizeof(double));
          for(unsigned int et=1; et<en; et++){ // loop over intermediate tau energies
            double dEt = delE[et-1];
            tau_decay_fluxes[et]    +=    nu_tau_flux_invlen_CC*    dNdE_CC_tau[et]*dEt;
            tau_bar_decay_fluxes[et]+=nu_tau_bar_flux_invlen_CC*dNdE_CC_tau_bar[et]*dEt;
          }
        }
      }
      
      //then accumulate the contributions back to the neutrino fluxes from the taus decaying
      
      //zero out the arrays
      std::fill(tau_hadlep_decays.begin(),tau_hadlep_decays.end(),0.);
      std::fill(tau_lep_decays.begin(),tau_lep_decays.end(),0.);
      //fill in new data
      for(unsigned int et=1; et<ne; et++){ // loop over intermediate tau energies
        double* tau_all_ptr=&int_struct->dNdE_tau_all[0][et][0];
        double* tau_lep_ptr=&int_struct->dNdE_tau_lep[0][et][0];
        double* taubar_all_ptr=&int_struct->dNdE_tau_all[1][et][0];
        double* taubar_lep_ptr=&int_struct->dNdE_tau_lep[1][et][0];
        SQUIDS_POINTER_IS_ALIGNED(tau_all_ptr,preferred_alignment*sizeof(double));
        SQUIDS_POINTER_IS_ALIGNED(tau_lep_ptr,preferred_alignment*sizeof(double));
        SQUIDS_POINTER_IS_ALIGNED(taubar_all_ptr,preferred_alignment*sizeof(double));
        SQUIDS_POINTER_IS_ALIGNED(taubar_lep_ptr,preferred_alignment*sizeof(double));
        double     tau_decay_flux_et=    tau_decay_fluxes[et];
        double tau_bar_decay_flux_et=tau_bar_decay_fluxes[et];
        for(unsigned int e1=0; e1<et; e1++){ // loop over final neutrino energies
          double tau_all_e1=tau_all_ptr[e1];
          double tau_lep_e1=tau_lep_ptr[e1];
          double taubar_all_e1=taubar_all_ptr[e1];
          double taubar_lep_e1=taubar_lep_ptr[e1];
          tau_hadlep_decays[0][e1] +=     tau_decay_flux_et*tau_all_e1;
          tau_lep_decays   [0][e1] += tau_bar_decay_flux_et*taubar_lep_e1;
          tau_hadlep_decays[1][e1] += tau_bar_decay_flux_et*taubar_all_e1;
          tau_lep_decays   [1][e1] +=     tau_decay_flux_et*tau_lep_e1;
        }
      }
      
    } //end of TR handling
    
    NUSQUIDS_DEBUG_LABEL(GR_contribution_osc);
    if(iglashow && (NT == both or NT == antineutrino)){
      std::fill(gr_factors.begin(),gr_factors.end(),0.);
      unsigned int rho=(NT == both) ? 1 : 0;
      
      for(unsigned int e2=1; e2<ne; e2++){
        double flux=evol_b1_proj[rho][0][e2]*estate[e2].rho[rho];
        double flux_invlen_en=flux*int_state.invlen_GR[e2]*delE[e2-1];
        double* dNdE_GR_ptr=&int_struct->dNdE_GR[e2][0];
        SQUIDS_POINTER_IS_ALIGNED(dNdE_GR_ptr,preferred_alignment*sizeof(double));
        for(unsigned int e1=0; e1<e2; e1++)
          gr_factors[e1] += flux_invlen_en*dNdE_GR_ptr[e1];
      }
    }
    
  }
  
  #undef ALIGNED_LOCAL_BUFFER
}

void nuSQUIDS::InitializeInteractions(){
  if(interactions_initialized)
    return; //nothing to do
  //===============================
  // init XS and TDecay objects  //
  //===============================

  try{
    // initialize tau decay spectra object
    tdc.Init(E_range);
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    throw std::runtime_error("nuSQUIDS::init : Failed while trying initialize TauDecaySpectra object [TauDecaySpectra::Init]");
  }
  
  // if we don't already have cached data, create it
  if ( !int_struct ) {
    int_struct.reset(new InteractionStructure);
    // initialize cross section object if we don't have one
    if ( ncs == nullptr) {
      ncs = std::make_shared<CrossSectionLibrary>(loadDefaultCrossSections());
    } // else we assume the user has already inintialized the object if not throw error.
    
    //We have to figure out what cross sections we are using. 
    //We would prefer to use separate ones for protons and neutrons, but it is 
    //possible that only isoscalar nuclear cross sections are available, in which
    //case we must use those. 
    //Note: This would be the starting point for supporting more advanced targets,
    //      e.g. whole isotopes, etc.
    bool perNucleonXS = ncs->hasTarget(proton) and ncs->hasTarget(neutron);
    unsigned int nTargets=0;
    if(perNucleonXS){
      int_struct->targets={proton,neutron};
      nTargets=2;
    }
    else{
      if(ncs->hasTarget(isoscalar_nucleon)){
        int_struct->targets={isoscalar_nucleon};
        nTargets=1;
      }
      else
        throw std::runtime_error("Cross section object does not provide any suitable nucleon cross sections");
    }
    
    // initialize cross section and interaction arrays
    try {
      InitializeInteractionVectors(nTargets);
    } catch (std::exception& ex) {
      std::cerr << ex.what() << std::endl;
      throw std::runtime_error("nuSQUIDS::init : Failed while trying to initialize interaction vectors [InitializeInteractionVectors]");
    }
    //===============================
    // Fill in arrays              //
    //===============================
    try {
      GetCrossSections();
    } catch (std::exception& ex) {
      std::cerr << ex.what() << std::endl;
      throw std::runtime_error("nuSQUIDS::init : Failed while trying to fill in interaction vectors [InitializeInteractions]");
    }
  }
  
  if(iinteraction){
    size_t rounded_ne=round_up_to_aligned(ne);
    nc_factors.resize(std::vector<size_t>{nrhos,3,rounded_ne});
    if(tauregeneration){
      tau_hadlep_decays.resize(std::vector<size_t>{2,rounded_ne});
      tau_lep_decays.resize(std::vector<size_t>{2,rounded_ne});
    }
    if(iglashow)
      gr_factors.resize(std::vector<size_t>{rounded_ne});
  }
  
  interactions_initialized=true;
}
  
void nuSQUIDS::GetCrossSections(){

    //units
    double cm2GeV = pow(params.cm,2)*pow(params.GeV,-1);
    double cm2 = pow(params.cm,2);
    double GeVm1 = pow(params.GeV,-1);

    // load cross sections

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

    auto validateCrossSection=[this](double value, double unit, const char* interaction, bool singleDiff, double neuE, double leptE, int flavor){
      //all active neutrino types should have sensible cross sections
      if(flavor<3 && (value<0.0 || std::isinf(value) || std::isnan(value))){
        std::ostringstream ss;
        ss << "Invalid " << interaction << (singleDiff ? " singly-differential" : " total")
        << " cross section value: " << value/unit << '\n';
        ss << " for flavor " << flavor << ",\n";
        ss << " neutrino energy " << neuE/this->params.GeV << " GeV";
        if(singleDiff)
          ss << ",\n outgoing lepton energy " << leptE/this->params.GeV << " GeV";
        throw std::runtime_error(ss.str());
      }
    };
    
//     std::cout << "Looking up cross sections. . . ";
//     std::cout.flush();
//     std::chrono::high_resolution_clock::time_point t1, t2;
//     t1 = std::chrono::high_resolution_clock::now();
  for(unsigned int target = 0; target < int_struct->targets.size(); target++){
    auto xs = ncs->crossSectionForTarget(int_struct->targets[target]);
    // initializing cross section arrays temporary array
    marray<double,4> dsignudE_CC{nrhos,numneu,ne,ne};
    marray<double,4> dsignudE_NC{nrhos,numneu,ne,ne};
    try{
    
    for(unsigned int neutype = 0; neutype < nrhos; neutype++){
      for(unsigned int flv = 0; flv < numneu; flv++){
        for(unsigned int e1 = 0; e1 < ne; e1++){
          // differential cross sections
          dsignudE_NC[neutype][flv][e1][0] = xs->SingleDifferentialCrossSection(E_range[e1],E_range[0],(NeutrinoCrossSections::NeutrinoFlavor)flv,neutype_xs_dict[neutype],NeutrinoCrossSections::NC)*cm2GeV;
          validateCrossSection(dsignudE_NC[neutype][flv][e1][0],cm2GeV,"NC",true,E_range[e1],E_range[0],flv);
          dsignudE_CC[neutype][flv][e1][0] = xs->SingleDifferentialCrossSection(E_range[e1],E_range[0],(NeutrinoCrossSections::NeutrinoFlavor)flv,neutype_xs_dict[neutype],NeutrinoCrossSections::CC)*cm2GeV;
          validateCrossSection(dsignudE_CC[neutype][flv][e1][0],cm2GeV,"CC",true,E_range[e1],E_range[0],flv);
          for(unsigned int e2 = 1; e2 < e1; e2++){
            dsignudE_NC[neutype][flv][e1][e2] = xs->AverageSingleDifferentialCrossSection(E_range[e1],E_range[e2-1],E_range[e2],(NeutrinoCrossSections::NeutrinoFlavor)flv,neutype_xs_dict[neutype],NeutrinoCrossSections::NC)*cm2GeV;
            validateCrossSection(dsignudE_NC[neutype][flv][e1][e2],cm2GeV,"NC",true,E_range[e1],E_range[e2],flv);
            dsignudE_CC[neutype][flv][e1][e2] = xs->AverageSingleDifferentialCrossSection(E_range[e1],E_range[e2-1],E_range[e2],(NeutrinoCrossSections::NeutrinoFlavor)flv,neutype_xs_dict[neutype],NeutrinoCrossSections::CC)*cm2GeV;
            validateCrossSection(dsignudE_CC[neutype][flv][e1][e2],cm2GeV,"CC",true,E_range[e1],E_range[e2],flv);
          }
          // total cross sections
          if(e1>0) {
            int_struct->sigma_CC[target][neutype][flv][e1] = xs->AverageTotalCrossSection(E_range[e1-1],E_range[e1],(NeutrinoCrossSections::NeutrinoFlavor)flv,neutype_xs_dict[neutype],NeutrinoCrossSections::CC)*cm2;
            int_struct->sigma_NC[target][neutype][flv][e1] = xs->AverageTotalCrossSection(E_range[e1-1],E_range[e1],(NeutrinoCrossSections::NeutrinoFlavor)flv,neutype_xs_dict[neutype],NeutrinoCrossSections::NC)*cm2;
          }
          else {
            int_struct->sigma_CC[target][neutype][flv][e1] = xs->TotalCrossSection(E_range[e1],static_cast<NeutrinoCrossSections::NeutrinoFlavor>(flv),neutype_xs_dict[neutype],NeutrinoCrossSections::CC)*cm2;
            int_struct->sigma_NC[target][neutype][flv][e1] = xs->TotalCrossSection(E_range[e1],static_cast<NeutrinoCrossSections::NeutrinoFlavor>(flv),neutype_xs_dict[neutype],NeutrinoCrossSections::NC)*cm2;
          }
          validateCrossSection(int_struct->sigma_CC[target][neutype][flv][e1],cm2,"CC",false,E_range[e1],0,flv);
          validateCrossSection(int_struct->sigma_NC[target][neutype][flv][e1],cm2,"NC",false,E_range[e1],0,flv);
        }
      }
    }

    // constructing dNdE for DIS
    for(unsigned int rho = 0; rho < nrhos; rho++){
      for(unsigned int flv = 0; flv < numneu; flv++){
        for(unsigned int e1 = 0; e1 < ne; e1++){
          for(unsigned int e2 = 0; e2 < e1; e2++){
            if(dsignudE_NC[rho][flv][e1][e2] == 0 )
              int_struct->dNdE_NC[target][rho][flv][e1][e2] = 0;
            else {
              int_struct->dNdE_NC[target][rho][flv][e1][e2] = (dsignudE_NC[rho][flv][e1][e2])/(int_struct->sigma_NC[target][rho][flv][e1]);
              validateCrossSection(int_struct->dNdE_NC[target][rho][flv][e1][e2],1.,"dNdE_NC",true,E_range[e1],E_range[e2],flv);
            }
            if(dsignudE_CC[rho][flv][e1][e2] == 0 )
              int_struct->dNdE_CC[target][rho][flv][e1][e2] = 0;
            else {
              int_struct->dNdE_CC[target][rho][flv][e1][e2] = (dsignudE_CC[rho][flv][e1][e2])/(int_struct->sigma_CC[target][rho][flv][e1]);
              validateCrossSection(int_struct->dNdE_CC[target][rho][flv][e1][e2],1.,"dNdE_CC",true,E_range[e1],E_range[e2],flv);
            }
          }
        }
      }
    }
      
    }catch(std::runtime_error& err){
      throw std::runtime_error(err.what()+std::string(" on target ")+std::to_string(int_struct->targets[target]));
    }
  }

    // construct dNdE for Glashow resonance
    if(iglashow && (NT == both or NT == antineutrino)){
      auto gr_cs = ncs->crossSectionForTarget(electron);
      if(!gr_cs)
        throw std::runtime_error("Glashow Resonance handling is active, but the supplied cross sections do not include a cross section for electrons");
      marray<double,2> dsignudE_GR{ne,ne};
      for(unsigned int e1 = 0; e1 < ne; e1++){
        if(e1>0)
          int_struct->sigma_GR[e1] = gr_cs->AverageTotalCrossSection(E_range[e1-1],E_range[e1],NeutrinoCrossSections::electron,NeutrinoCrossSections::antineutrino,NeutrinoCrossSections::GR)*cm2;
        else
          int_struct->sigma_GR[e1] = gr_cs->TotalCrossSection(E_range[e1],NeutrinoCrossSections::electron,NeutrinoCrossSections::antineutrino,NeutrinoCrossSections::GR)*cm2;
        dsignudE_GR[e1][0] = gr_cs->SingleDifferentialCrossSection(E_range[e1],E_range[0],NeutrinoCrossSections::electron,NeutrinoCrossSections::antineutrino,NeutrinoCrossSections::GR)*cm2GeV;
        for(unsigned int e2 = 1; e2 < e1; e2++)
          dsignudE_GR[e1][e2] = gr_cs->AverageSingleDifferentialCrossSection(E_range[e1],E_range[e2-1],E_range[e2],NeutrinoCrossSections::electron,NeutrinoCrossSections::antineutrino,NeutrinoCrossSections::GR)*cm2GeV;
      }
      for(unsigned int e1 = 0; e1 < ne; e1++){
        for(unsigned int e2 = 0; e2 < e1; e2++){
          int_struct->dNdE_GR[e1][e2] = dsignudE_GR[e1][e2]/int_struct->sigma_GR[e1];
        }
      }
    }

    // load tau decay spectra
    //
    // filling cross section arrays
    std::map<unsigned int,unsigned int> neutype_decay_dict;
    if (NT == neutrino){
      neutype_decay_dict = (std::map<unsigned int,unsigned int>){{0,0}};
    } else if ( NT == antineutrino ) {
      neutype_decay_dict = (std::map<unsigned int,unsigned int>){{0,1}};
    } else {
      // in this case NT is both
      neutype_decay_dict = (std::map<unsigned int,unsigned int>){{0,0},{1,1}};
    }
    // constructing dNdE_tau_lep/dNdE_tau_all
    for(unsigned int neutype = 0; neutype < nrhos; neutype++){
      for(unsigned int e1 = 0; e1 < ne; e1++){
          for(unsigned int e2 = 0; e2 < e1; e2++){
              int_struct->dNdE_tau_all[neutype][e1][e2] = tdc.dNdEnu_All(e1,e2,neutype_decay_dict[neutype])*GeVm1;
              int_struct->dNdE_tau_lep[neutype][e1][e2] = tdc.dNdEnu_Lep(e1,e2,neutype_decay_dict[neutype])*GeVm1;
          }
      }
    }
    
//     t2 = std::chrono::high_resolution_clock::now();
//     double time=std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
//     std::cout << time << " seconds" << std::endl;
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
          estate[ie].rho[rho] -= evol_b1_proj[rho][flv][ie]*quantity;
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
  if ( body == nullptr )
    throw std::runtime_error("nuSQUIDS::Error::BODY is a NULL pointer");
  if (not ibody )
    throw std::runtime_error("nuSQUIDS::Error::Body not initialized");
  if ( track == nullptr )
    throw std::runtime_error("nuSQUIDS::Error::TRACK is a NULL pointer");
  if ( not itrack )
    throw std::runtime_error("nuSQUIDS::Error::TRACK is not initialized");
  if ( not istate )
    throw std::runtime_error("nuSQUIDS::Error::Initial state not initialized");
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  
  if( iinteraction && !interactions_initialized )
    InitializeInteractions();
  if( !ioscillations && iinteraction)
    SetUpInteractionCache();

  // is track time reversed
  if(track->GetFinalX() < track->GetInitialX()){
    // flip the arrow of time
    Set_h((-1.0)*Get_h());
  }

  if ( (body->IsConstantDensity() and not iinteraction) and allowConstantDensityOscillationOnlyEvolution ){
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
        //here we deliberately update state directly,
        //rather than working on estate as we would is using GSL ODE evolution
        tmp2 = state[ie].rho[rho].UTransform(tmp1,gsl_complex_rect(0.,evolution_time));
        state[ie].rho[rho] = tmp2;
        //std::cout << "despues " << E_range[ie] << " " << tmp2 << std::endl;
      }
    }
    if(progressbar)
      ProgressBar();
    return;
  }

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
  if(progressbar)
    ProgressBar();
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
      state[i].rho[r].SetAllComponents(0);
      if (basis == flavor){
        for(unsigned int j = 0; j < v.extent(0); j++)
          state[i].rho[r] += v[j]*b1_proj[r][j];
      }
      else if (basis == mass){
        for(int j = 0; j < v.extent(0); j++)
          state[i].rho[r] += v[j]*b0_proj[j];
      }
    }
  }

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
      state[i].rho[r].SetAllComponents(0);
      if (basis == flavor){
        for(unsigned int j = 0; j < numneu; j++)
          state[i].rho[r] += v[i][j]*b1_proj[r][j];
      }
      else if (basis == mass){
        for(unsigned int j = 0; j < numneu; j++)
          state[i].rho[r] += v[i][j]*b0_proj[j];
      }
    }
  }

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
      state[i].rho[r].SetAllComponents(0);
      if (basis == flavor){
        for(unsigned int j = 0; j < numneu; j++)
          state[i].rho[r] += v[i][r][j]*b1_proj[r][j];
      }
      else if (basis == mass){
        for(unsigned int j = 0; j < numneu; j++)
          state[i].rho[r] += v[i][r][j]*b0_proj[j];
      }
    }
  }
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
    return b0_proj[flv]*GetIntermediateState(rho,EE);
  if ( EE < *E_range.begin() || EE > *E_range.rbegin() )
    throw std::runtime_error("nuSQUIDS::Error::Energy "+std::to_string(EE)+" outside of propagated energy range, ["
                             +std::to_string(*E_range.begin())+","+std::to_string(*E_range.rbegin())+"].");
  return GetExpectationValueD(b0_proj[flv], rho, EE);
}

double nuSQUIDS::EvalFlavor(unsigned int flv,double EE,unsigned int rho) const{
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != both )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  if ( EE < *E_range.begin() || EE > *E_range.rbegin() )
    throw std::runtime_error("nuSQUIDS::Error::Energy "+std::to_string(EE)+" outside of propagated energy range, ["
                             +std::to_string(*E_range.begin())+","+std::to_string(*E_range.rbegin())+"].");
  if ( basis == mass )
    return b1_proj[rho][flv]*GetIntermediateState(rho,EE);
  return GetExpectationValueD(b1_proj[rho][flv], rho, EE);
}

double nuSQUIDS::EvalMass(unsigned int flv,double EE, unsigned int rho, double scale, std::vector<bool>& avr) const{
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != both )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  if ( EE < *E_range.begin() || EE > *E_range.rbegin() )
    throw std::runtime_error("nuSQUIDS::Error::Energy "+std::to_string(EE)+" outside of propagated energy range, ["
                             +std::to_string(*E_range.begin())+","+std::to_string(*E_range.rbegin())+"].");
  if (avr.size() != numneu)
    throw std::runtime_error("nuSQUIDS::Error::Average vector bool size must be the number of neutrinos.");
  if ( basis == mass )
    return b0_proj[flv]*GetIntermediateState(rho,EE);
  return GetExpectationValueD(b0_proj[flv], rho, EE, scale, avr);
}

double nuSQUIDS::EvalFlavor(unsigned int flv,double EE,unsigned int rho, double scale, std::vector<bool>& avr) const{
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");
  if ( rho != 0 and NT != both )
    throw std::runtime_error("nuSQUIDS::Error::Cannot evaluate rho != 0 in this NT mode.");
  if ( EE < *E_range.begin() || EE > *E_range.rbegin() )
    throw std::runtime_error("nuSQUIDS::Error::Energy "+std::to_string(EE)+" outside of propagated energy range, ["
                             +std::to_string(*E_range.begin())+","+std::to_string(*E_range.rbegin())+"].");
  if (avr.size() != numneu)
    throw std::runtime_error("nuSQUIDS::Error::Average vector bool size must be the number of neutrinos.");
  if ( basis == mass )
    return b1_proj[rho][flv]*GetIntermediateState(rho,EE);
  return GetExpectationValueD(b1_proj[rho][flv], rho, EE, scale, avr);
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

double nuSQUIDS::EvalMass(unsigned int flv, double scale, std::vector<bool>& avr) const{
  if(state == nullptr)
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
  return GetExpectationValue(b0_proj[flv], 0, 0, scale, avr);
}

double nuSQUIDS::EvalFlavor(unsigned int flv, double scale, std::vector<bool>& avr) const{
  if(state == nullptr)
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
  return GetExpectationValue(b1_proj[0][flv], 0, 0, scale, avr);
}

double nuSQUIDS::EvalMass(unsigned int flv) const{
  if(state == nullptr)
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
  if(state == nullptr)
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

squids::SU_vector nuSQUIDS::GetHamiltonian(unsigned int ei, unsigned int rho){
  if (!ienergy)
    throw std::runtime_error("nuSQUIDS::Error::Energy not initialized");
  PreDerive(Get_t());
  return H0(E_range[ei],rho)+HI(ei,rho,Get_t());
}

void nuSQUIDS::WriteStateHDF5(std::string str,std::string grp,bool save_cross_section, std::string cross_section_grp_loc, bool overwrite) const{
  if ( body == nullptr )
    throw std::runtime_error("nuSQUIDS::Error::BODY is a NULL pointer");
  if (not ibody )
    throw std::runtime_error("nuSQUIDS::Error::Body not initialized");
  if ( track == nullptr )
    throw std::runtime_error("nuSQUIDS::Error::TRACK is a NULL pointer");
  if ( not itrack )
    throw std::runtime_error("nuSQUIDS::Error::TRACK is not initialized");
  if ( not istate )
    throw std::runtime_error("nuSQUIDS::Error::Initial state not initialized");
  if ( not ienergy )
    throw std::runtime_error("nuSQUIDS::Error::Energy not set.");

  if (!iinteraction)
    save_cross_section = false;

  // this lines supress HDF5 error messages
  H5Eset_auto (H5E_DEFAULT,nullptr, nullptr);

  //hid_t dset_id;
  // create HDF5 file
  // H5F_ACC_TRUNC : overwrittes file
  // H5F_ACC_EXCL  : files if file exists
  H5File file(overwrite ?
              H5Fcreate(str.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT) :
              H5Fopen(str.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));

  if (file < 0) {// file could not be open. Create a new file.
    file=H5File(H5Fcreate(str.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    if (file < 0)
        throw std::runtime_error("nuSQUIDS: Error: Cannot create file at " + str + ".");
  }

  H5Handle rootGroup(H5Gopen(file, "/",H5P_DEFAULT), H5Gclose, "open HDF5 group");
  if (strncmp(grp.c_str(),"/",1)!=0){
    std::cout << "nuSQUIDS::WriteStateHDF5::Warning::group location name did not start with '/'. '/' was prepend to " + grp << std::endl;
    grp = "/"+grp;
  }

  H5Handle group;
  if ( grp != "/" )
    group = H5Handle(H5Gcreate(rootGroup, grp.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),H5Gclose,"create HDF5 group");
  else
    group = H5Handle(rootGroup.get(),[](hid_t)->herr_t{ return 0; },"get HDF5 group"); //non-owning handle

  // write the energy range
  hsize_t Edims[1]={E_range.extent(0)};
  H5LTmake_dataset(group,"energies",1,Edims,H5T_NATIVE_DOUBLE,E_range.get_data());

  // write mixing parameters
  hsize_t dim[1]{1};
  H5LTmake_dataset(group,"basic",1,dim,H5T_NATIVE_DOUBLE,nullptr);
  H5LTmake_dataset(group,"mixingangles",1,dim,H5T_NATIVE_DOUBLE,nullptr);
  H5LTmake_dataset(group,"CPphases",1,dim,H5T_NATIVE_DOUBLE,nullptr);
  H5LTmake_dataset(group,"massdifferences",1,dim,H5T_NATIVE_DOUBLE,nullptr);

  H5LTset_attribute_int(group, "basic","numneu",(const int*)&numneu, 1);
  int auxint = static_cast<int>(NT);
  H5LTset_attribute_int(group, "basic","NT",&auxint,1);
  H5LTset_attribute_string(group, "basic", "interactions", (iinteraction) ? "True":"False");
  auxint = static_cast<int>(ioscillations);
  H5LTset_attribute_int(group, "basic", "oscillations", &auxint, 1);
  auxint = static_cast<int>(tauregeneration);
  H5LTset_attribute_int(group, "basic", "tau_regeneration", &auxint, 1);
  auxint = static_cast<int>(iglashow);
  H5LTset_attribute_int(group, "basic", "glashow_resonance", &auxint, 1);
  double auxt = Get_t();
  H5LTset_attribute_double(group, "basic", "squids_time", &auxt,1);
  double auxt_ini = Get_t_initial();
  H5LTset_attribute_double(group, "basic", "squids_time_initial", &auxt_ini,1);

  // version numbers
  H5LTset_attribute_string(group, "basic", "squids_version", SQUIDS_VERSION_STR);
  unsigned int squids_version = SQUIDS_VERSION;
  H5LTset_attribute_uint(group, "basic", "squids_version_number", &squids_version,1);

  H5LTset_attribute_string(group, "basic", "nusquids_version", NUSQUIDS_VERSION_STR);
  unsigned int nusquids_version = NUSQUIDS_VERSION;
  H5LTset_attribute_uint(group, "basic", "nusquids_version_number", &nusquids_version,1);

  // set mixing angles
  for( unsigned int i = 0; i < numneu; i++ ){
    for( unsigned int j = i+1; j < numneu; j++ ){
      std::string th_label = "th"+std::to_string(i+1)+std::to_string(j+1);
      double th_value = params.GetMixingAngle(i,j);
      H5LTset_attribute_double(group, "mixingangles",th_label.c_str(),&th_value, 1);

      std::string delta_label = "delta"+std::to_string(i+1)+std::to_string(j+1);
      double delta_value = params.GetPhase(i,j);
      H5LTset_attribute_double(group, "CPphases",delta_label.c_str(),&delta_value, 1);
    }
  }

  for ( unsigned int i = 1; i < numneu; i++ ){
    std::string dm2_label = "dm"+std::to_string(i+1)+"1sq";
    double dm2_value = params.GetEnergyDifference(i);
    H5LTset_attribute_double(group, "massdifferences",dm2_label.c_str(),&dm2_value, 1);
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

  H5LTmake_dataset(group,"neustate",2,statedim,H5T_NATIVE_DOUBLE,static_cast<const void*>(neustate.data()));
  H5LTmake_dataset(group,"aneustate",2,statedim,H5T_NATIVE_DOUBLE,static_cast<const void*>(aneustate.data()));

  // writing body and track information
  hid_t track_group_id = H5Gcreate(group, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track->Serialize(track_group_id);
  H5Gclose(track_group_id);

  hid_t body_group_id = H5Gcreate(group, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body->Serialize(body_group_id);
  H5Gclose(body_group_id);

  // writing cross section information
  H5Handle xs_group;
  if ( cross_section_grp_loc == ""){
    xs_group = H5Handle(H5Gcreate(group, "crosssections", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),H5Gclose," create HDF5 group");
  } else {
    //this location is shared, so the group may already exist
    if(h5ObjectExists(rootGroup, cross_section_grp_loc.c_str())){
      xs_group = H5Handle(H5Gopen(rootGroup, cross_section_grp_loc.c_str(), H5P_DEFAULT), H5Gclose, "open HDF5 shared crosssection group");
    }
    else{ //but if it does not exist, just create it
      xs_group = H5Handle(H5Gcreate(rootGroup, cross_section_grp_loc.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),H5Gclose,"create HDF5 shared crosssection group");
    }
  }

  if(iinteraction and save_cross_section) {
    if(!interactions_initialized){
      // we shall pay for this sin. Evil. Very evil.
      const_cast<nuSQUIDS*>(this)->InitializeInteractions();
    }
    
    auto nTargets = int_struct->targets.size();
    { //target type list
      static_assert(sizeof(std::underlying_type<PDGCode>::type)==(32/CHAR_BIT), "Following code assumes PDGCode is int32_t");
      static_assert(std::is_signed<std::underlying_type<PDGCode>::type>::value, "Following code assumes PDGCode is int32_t");
      
      hsize_t rawNTargets=nTargets;
      H5Handle dspaceTargets(H5Screate_simple(1,&rawNTargets,nullptr), H5Sclose, "create HDF5 dataspace");
      H5Handle dsetTargets(H5Dcreate2(xs_group,"targets",H5T_NATIVE_INT32,dspaceTargets,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT), H5Dclose, "create HDF5 dataset");
      
      herr_t err=H5Dwrite(dsetTargets, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->targets.begin());
      if(err<0)
        throw std::runtime_error("Failed to write data to HDF5 dataset (targets)");
      addH5Attribute<std::string>(dsetTargets,"Description","PDG Codes for target matter types");
    }
    { // total crosssections
      hsize_t XSdim[4] {nTargets, nrhos, numneu, ne};
        
      H5Handle dspaceSigma(H5Screate_simple(4,XSdim,nullptr), H5Sclose, "create HDF5 dataspace");
      H5Handle dsetSigmaCC(H5Dcreate2(xs_group,"sigmacc",H5T_NATIVE_DOUBLE,dspaceSigma,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT), H5Dclose, "create HDF5 dataset");
      H5Handle dsetSigmaNC(H5Dcreate2(xs_group,"sigmanc",H5T_NATIVE_DOUBLE,dspaceSigma,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT), H5Dclose, "create HDF5 dataset");
      H5Handle dspaceSigmaGR(H5Screate_simple(1,XSdim+3,nullptr), H5Sclose, "create HDF5 dataspace"); //GR has lower rank
      H5Handle dsetSigmaGR(H5Dcreate2(xs_group,"sigmagr",H5T_NATIVE_DOUBLE,dspaceSigmaGR,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT), H5Dclose, "create HDF5 dataset");
      
      //check the last dimension of the cross section arrays for padding
      if(int_struct->sigma_CC.extent(3)==ne){
        //no padding, which makes writing easy; dataspaces are the same shape
        //write sigma CC
        herr_t err=H5Dwrite(dsetSigmaCC, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->sigma_CC.begin());
        if(err<0)
          throw std::runtime_error("Failed to write data to HDF5 dataset (sigma CC)");
        
        //write sigma NC
        err=H5Dwrite(dsetSigmaNC, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->sigma_NC.begin());
        if(err<0)
          throw std::runtime_error("Failed to write data to HDF5 dataset (sigma NC)");
        
        //write sigma NC
        err=H5Dwrite(dsetSigmaGR, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->sigma_GR.begin());
        if(err<0)
          throw std::runtime_error("Failed to write data to HDF5 dataset (sigma GR)");
      }
      else{
        //after each run of data in the last dimension, the array has padding entries for alignment
        //which we do not want to write. Instead, we write the blocks which fit between the 
        //sections of padding. 
        
        //The in-memory dataspace is a different size. 
        hsize_t source_dim[4] {int_struct->sigma_CC.extent(0),
          int_struct->sigma_CC.extent(1),
          int_struct->sigma_CC.extent(2),
          int_struct->sigma_CC.extent(3)};
        H5Handle mem_dspace(H5Screate_simple(4,source_dim,nullptr), H5Sclose, "create HDF5 dataspace");
        H5Handle mem_dspaceGR(H5Screate_simple(1,source_dim+3,nullptr), H5Sclose, "create HDF5 dataspace");
        
        hsize_t stride[4]={1,1,1,1}; //we want everything, so all strides are minimal
        hsize_t count[4]={nTargets,nrhos,numneu,1}; //step through all but the last dimension taking blocks
        hsize_t block[4]={1,1,1,ne}; //the individual blocks are this size
        hsize_t start[4]={0,0,0,0};
        //designate which section of the memory we will read, and the file contents we will write
        herr_t err=H5Sselect_hyperslab(dspaceSigma, H5S_SELECT_SET, start, stride, count, block);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 on-disk dataspace");
        err=H5Sselect_hyperslab(mem_dspace.get(), H5S_SELECT_SET, start, stride, count, block);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 in-memory dataspace");
        err=H5Sselect_hyperslab(dspaceSigmaGR, H5S_SELECT_SET, start+3, stride+3, count+3, block+3);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 on-disk dataspace");
        err=H5Sselect_hyperslab(mem_dspaceGR.get(), H5S_SELECT_SET, start+3, stride+3, count+3, block+3);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 in-memory dataspace");
        
        //write sigma CC
        err=H5Dwrite(dsetSigmaCC, H5T_NATIVE_DOUBLE, mem_dspace, dspaceSigma, H5P_DEFAULT, &*int_struct->sigma_CC.begin());
        if(err<0)
          throw std::runtime_error("Failed to write masked data to HDF5 dataset (sigma CC)");
        
        //write sigma NC
        err=H5Dwrite(dsetSigmaNC, H5T_NATIVE_DOUBLE, mem_dspace, dspaceSigma, H5P_DEFAULT, &*int_struct->sigma_NC.begin());
        if(err<0)
          throw std::runtime_error("Failed to write masked data to HDF5 dataset (sigma NC)");
        
        //write sigma GR
        err=H5Dwrite(dsetSigmaGR, H5T_NATIVE_DOUBLE, mem_dspaceGR, dspaceSigmaGR, H5P_DEFAULT, &*int_struct->sigma_GR.begin());
        if(err<0)
          throw std::runtime_error("Failed to write masked data to HDF5 dataset (sigma GR)");
      }
      
      addH5Attribute<std::string>(dsetSigmaCC,"Description","Tabulated total cross sections for charged-current interactions");
      addH5Attribute<std::string>(dsetSigmaNC,"Description","Tabulated total cross sections for neutral-current interactions");
      addH5Attribute<std::string>(dsetSigmaGR,"Description","Tabulated total cross sections for Glashow resonance interactions");
    }

    { // differential crossections
      hsize_t dXSdim[5] {nTargets, nrhos, numneu, ne, ne};
      
      H5Handle dspaceDSigma(H5Screate_simple(5,dXSdim,nullptr), H5Sclose, "create HDF5 dataspace");
      H5Handle dsetDSigmaCC(H5Dcreate2(xs_group,"dNdEcc",H5T_NATIVE_DOUBLE,dspaceDSigma,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT), H5Dclose, "create HDF5 dataset");
      H5Handle dsetDSigmaNC(H5Dcreate2(xs_group,"dNdEnc",H5T_NATIVE_DOUBLE,dspaceDSigma,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT), H5Dclose, "create HDF5 dataset");
      H5Handle dspaceDSigmaGR(H5Screate_simple(2,dXSdim+3,nullptr), H5Sclose, "create HDF5 dataspace"); //GR has lower rank
      H5Handle dsetDSigmaGR(H5Dcreate2(xs_group,"dNdEgr",H5T_NATIVE_DOUBLE,dspaceDSigmaGR,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT), H5Dclose, "create HDF5 dataset");
      
      //check the last dimension of the cross section arrays for padding
      if(int_struct->dNdE_GR.extent(4)==ne){
        //no padding, which makes writing easy; dataspaces are the same shape
        //write dNdE CC
        herr_t err=H5Dwrite(dsetDSigmaCC, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_CC.begin());
        if(err<0)
          throw std::runtime_error("Failed to write data to HDF5 dataset (dsigma CC)");
        
        //write dNdE NC
        err=H5Dwrite(dsetDSigmaNC, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_NC.begin());
        if(err<0)
          throw std::runtime_error("Failed to write data to HDF5 dataset (dsigma NC)");
        
        //write dNdE NC
        err=H5Dwrite(dsetDSigmaGR, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_GR.begin());
        if(err<0)
          throw std::runtime_error("Failed to write data to HDF5 dataset (dsigma GR)");
      }
      else{ //harder case with padding
        //The in-memory dataspace is a different size. 
        hsize_t source_dim[5] {int_struct->dNdE_CC.extent(0),
          int_struct->dNdE_CC.extent(1),
          int_struct->dNdE_CC.extent(2),
          int_struct->dNdE_CC.extent(3),
          int_struct->dNdE_CC.extent(4)};
        H5Handle mem_dspace(H5Screate_simple(5,source_dim,nullptr), H5Sclose, "create HDF5 dataspace");
        H5Handle mem_dspaceGR(H5Screate_simple(2,source_dim+3,nullptr), H5Sclose, "create HDF5 dataspace");
        
        hsize_t stride[5]={1,1,1,1,1}; //we want everything, so all strides are minimal
        hsize_t count[5]={nTargets,nrhos,numneu,ne,1}; //step through all but the last dimension taking blocks
        hsize_t block[5]={1,1,1,1,ne}; //the individual blocks are this size
        hsize_t start[5]={0,0,0,0,0};
        //designate which section of the memory we will read, and the file contents we will write
        herr_t err=H5Sselect_hyperslab(dspaceDSigma, H5S_SELECT_SET, start, stride, count, block);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 on-disk dataspace");
        err=H5Sselect_hyperslab(mem_dspace, H5S_SELECT_SET, start, stride, count, block);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 in-memory dataspace");
        err=H5Sselect_hyperslab(dspaceDSigmaGR, H5S_SELECT_SET, start+3, stride+3, count+3, block+3);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 on-disk dataspace");
        err=H5Sselect_hyperslab(mem_dspaceGR, H5S_SELECT_SET, start+3, stride+3, count+3, block+3);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 in-memory dataspace");
        
        //write dNdE CC
        err=H5Dwrite(dsetDSigmaCC, H5T_NATIVE_DOUBLE, mem_dspace, dspaceDSigma, H5P_DEFAULT, &*int_struct->dNdE_CC.begin());
        if(err<0)
          throw std::runtime_error("Failed to write masked data to HDF5 dataset (dsigma CC)");
        
        //write dNdE NC
        err=H5Dwrite(dsetDSigmaNC, H5T_NATIVE_DOUBLE, mem_dspace, dspaceDSigma, H5P_DEFAULT, &*int_struct->dNdE_NC.begin());
        if(err<0)
          throw std::runtime_error("Failed to write masked data to HDF5 dataset (dsigma NC)");
        
        //write dNdE GR
        err=H5Dwrite(dsetDSigmaGR, H5T_NATIVE_DOUBLE, mem_dspaceGR, dspaceDSigmaGR, H5P_DEFAULT, &*int_struct->dNdE_GR.begin());
        if(err<0)
          throw std::runtime_error("Failed to write masked data to HDF5 dataset (dsigma GR)");
      }
      
      addH5Attribute<std::string>(dsetDSigmaCC,"Description","Tabulated differential cross sections for charged-current interactions");
      addH5Attribute<std::string>(dsetDSigmaNC,"Description","Tabulated differential cross sections for neutral-current interactions");
      addH5Attribute<std::string>(dsetDSigmaGR,"Description","Tabulated differential cross sections for Glashow resonance interactions");
    }

    { // tau decay data
      hsize_t dNdEtaudim[3] {nrhos, ne, ne};
      
      H5Handle dspaceTau(H5Screate_simple(3,dNdEtaudim,nullptr), H5Sclose, "create HDF5 dataspace");
      H5Handle dNdEtauall(H5Dcreate2(xs_group,"dNdEtauall",H5T_NATIVE_DOUBLE,dspaceTau,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT), H5Dclose, "create HDF5 dataset");
      H5Handle dNdEtaulep(H5Dcreate2(xs_group,"dNdEtaulep",H5T_NATIVE_DOUBLE,dspaceTau,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT), H5Dclose, "create HDF5 dataset");
      
      //check the last dimension of the cross section arrays for padding
      if(int_struct->dNdE_tau_all.extent(2)==ne){
        //no padding, which makes writing easy; dataspaces are the same shape
        herr_t err=H5Dwrite(dNdEtauall, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_tau_all.begin());
        if(err<0)
          throw std::runtime_error("Failed to write data to HDF5 dataset (dNdEtauall)");
        
        err=H5Dwrite(dNdEtaulep, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_tau_lep.begin());
        if(err<0)
          throw std::runtime_error("Failed to write data to HDF5 dataset (dNdEtaulep)");
      }
      else{ //harder case with padding
        //The in-memory dataspace is a different size. 
        hsize_t source_dim[3] {int_struct->dNdE_tau_all.extent(0),
          int_struct->dNdE_tau_all.extent(1),
          int_struct->dNdE_tau_all.extent(2)};
        H5Handle mem_dspace(H5Screate_simple(3,source_dim,nullptr), H5Sclose, "create HDF5 dataspace");
        
        hsize_t stride[3]={1,1,1}; //we want everything, so all strides are minimal
        hsize_t count[3]={nrhos,ne,1}; //step through all but the last dimension taking blocks
        hsize_t block[3]={1,1,ne}; //the individual blocks are this size
        hsize_t start[3]={0,0,0};
        //designate which section of the memory we will read, and the file contents we will write
        herr_t err=H5Sselect_hyperslab(dspaceTau, H5S_SELECT_SET, start, stride, count, block);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 on-disk dataspace");
        err=H5Sselect_hyperslab(mem_dspace, H5S_SELECT_SET, start, stride, count, block);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 in-memory dataspace");
        
        err=H5Dwrite(dNdEtauall, H5T_NATIVE_DOUBLE, mem_dspace, dspaceTau, H5P_DEFAULT, &*int_struct->dNdE_tau_all.begin());
        if(err<0)
          throw std::runtime_error("Failed to write masked data to HDF5 dataset (dNdEtauall)");
        
        err=H5Dwrite(dNdEtaulep, H5T_NATIVE_DOUBLE, mem_dspace, dspaceTau, H5P_DEFAULT, &*int_struct->dNdE_tau_lep.begin());
        if(err<0)
          throw std::runtime_error("Failed to write masked data to HDF5 dataset (dNdEtaulep)");
      }
      
      addH5Attribute<std::string>(dNdEtauall,"Description","Tabulated tau decay specrum for all channels");
      addH5Attribute<std::string>(dNdEtaulep,"Description","Tabulated tau decay specrum for leptonic channels");
    }
  }

  // write user parameters
  H5Handle user_parameters(H5Gcreate(group, "user_parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),H5Gclose,"create HDF5 group");
  // give control to the user and temporary restore HDF5 error messages
  //H5Eset_auto (H5E_DEFAULT,(H5E_auto_t) H5Eprint,stderr);
  AddToWriteHDF5(user_parameters);
  //H5Eset_auto (H5E_DEFAULT,nullptr,nullptr);
}

void nuSQUIDS::AddToWriteHDF5(hid_t hdf5_loc_id) const {

}

void nuSQUIDS::AddToReadHDF5(hid_t hdf5_loc_id){

}

void nuSQUIDS::ReadStateHDF5Internal(std::string str,std::string grp,std::shared_ptr<InteractionStructure> iis){
  // open HDF5 file
  //std::cout << "reading from hdf5 file" << std::endl;
  H5File file(H5Fopen(str.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (file < 0)
      throw std::runtime_error("nuSQUIDS::Error::file not found : " + str + ".");
  H5Handle rootGroup(H5Gopen(file, "/", H5P_DEFAULT), H5Gclose, "open HDF5 group");
  if (strncmp(grp.c_str(),"/",1)!=0){
    std::cout << "nuSQUIDS::ReadStateHDF5::Warning::group location name did not start with '/'. '/' was prepend to " + grp << std::endl;
    grp = "/"+grp;
  }
  H5Handle group(H5Gopen(rootGroup, grp.c_str(), H5P_DEFAULT),H5Gclose,"open HDF5 group");
  if ( group < 0 )
      throw std::runtime_error("nuSQUIDS::Error::Group '" + grp + "' does not exist in HDF5.");

  // read number of neutrinos
  H5LTget_attribute_uint(group, "basic", "numneu", static_cast<unsigned int*>(&numneu));
  // neutrino/antineutrino/both
  int auxint;
  char auxchar[20];
  herr_t err;
  H5LTget_attribute_int(group, "basic", "NT", &auxint);
  NT = static_cast<NeutrinoType>(auxint);
  // interactions
  H5LTget_attribute_string(group,"basic","interactions", auxchar);
  std::string aux = auxchar;
  if ( aux == "True")
    iinteraction = true;
  else
    iinteraction = false;

  err=H5LTget_attribute_int(group, "basic", "oscillations", &auxint);
  if(err>=0) ioscillations=auxint;
  err=H5LTget_attribute_int(group, "basic", "tau_regeneration", &auxint);
  if(err>=0) tauregeneration=auxint;
  err=H5LTget_attribute_int(group, "basic", "glashow_resonance", &auxint);
  if(err>=0) iglashow=auxint;

  double squids_time;
  H5LTget_attribute_double(group, "basic", "squids_time", &squids_time);

  double squids_time_initial;
  H5LTget_attribute_double(group, "basic", "squids_time_initial", &squids_time_initial);

  // check version numbers
  unsigned int squids_version;
  H5LTget_attribute_uint(group, "basic", "squids_version_number", &squids_version);

  if ( squids_version > SQUIDS_VERSION )
    throw std::runtime_error("nuSQUIDS::ReadStateHDF5::Error: File was written using SQuIDS version " +
        std::to_string(squids_version) + " current version is " + std::to_string(SQUIDS_VERSION));

  unsigned int nusquids_version;
  H5LTget_attribute_uint(group, "basic", "nusquids_version_number", &nusquids_version);
  if ( nusquids_version > NUSQUIDS_VERSION )
    throw std::runtime_error("nuSQUIDS::ReadStateHDF5::Error: File was written using nuSQuIDS version " +
        std::to_string(nusquids_version) + " current version is " + std::to_string(NUSQUIDS_VERSION));

  // reading body and track
  if(nusquids_version>100000){
    std::string body_name = getStringAttribute(group,"body","name");
    H5Handle body_group(H5Gopen(group, "body", H5P_DEFAULT), H5Gclose, "open body group");
    body=GetBodyDeserializer(body_name)(body_group);

    H5Handle track_group(H5Gopen(group, "track", H5P_DEFAULT), H5Gclose, "open track group");
    std::string track_name = getStringAttribute(group,"track","name");
    track=GetTrackDeserializer(track_name)(track_group);
  } else {
    unsigned int body_id;
    hsize_t dimbody[1];
    H5LTget_attribute_uint(group,"body","ID",&body_id);

    H5LTget_dataset_info(group,"body", dimbody,nullptr,nullptr);
    double body_params[dimbody[0]];
    H5LTread_dataset_double(group,"body", body_params);

    hsize_t dimtrack[1];
    H5LTget_dataset_info(group,"track", dimtrack ,nullptr,nullptr);
    std::unique_ptr<double[]> track_params(new double[dimtrack[0]]);
    H5LTread_dataset_double(group,"track", track_params.get());

    double x_current;
    H5LTget_attribute_double(group,"track","X",&x_current);

    // setting body and track
    SetBodyTrack(body_id,dimbody[0],body_params,dimtrack[0],track_params.get());
    // set trayectory to current time
    track->SetX(x_current);
  }

  // reading energy
  hsize_t dims[2];
  H5LTget_dataset_info(group, "energies", dims, nullptr, nullptr);

  ne = static_cast<unsigned int>(dims[0]);
  std::unique_ptr<double[]> energy_data(new double[ne]);
  H5LTread_dataset_double(group, "energies", energy_data.get());
  E_range = marray<double,1>{ne};
  for (unsigned int ie = 0; ie < ne; ie++)
    E_range[ie] = energy_data[ie];

  // initializing nuSQUIDS
  if (ne == 1){
    if(not inusquids)
      init(squids_time_initial);
    Set_E(energy_data[0]);
  }
  else {
    init(E_range,squids_time_initial);
  }


  // read and set mixing parameters
  for( unsigned int i = 0; i < numneu; i++ ){
    for( unsigned int j = i+1; j < numneu; j++ ){
      double th_value;
      std::string th_label = "th"+std::to_string(i+1)+std::to_string(j+1);
      H5LTget_attribute_double(group,"mixingangles", th_label.c_str(), &th_value);
      Set_MixingAngle(i,j,th_value);

      double delta_value;
      std::string delta_label = "delta"+std::to_string(i+1)+std::to_string(j+1);
      H5LTget_attribute_double(group,"CPphases", delta_label.c_str(), &delta_value);
      Set_CPPhase(i,j,delta_value);
    }
  }

  for( unsigned int i = 1; i < numneu; i++ ){
    double dm2_value;
    std::string dm2_label = "dm"+std::to_string(i+1)+"1sq";
    H5LTget_attribute_double(group,"massdifferences", dm2_label.c_str(), &dm2_value);
    Set_SquareMassDifference(i, dm2_value);
  }

  // reset current squids time
  Set_t(squids_time);
  // set time offset
  time_offset = squids_time - track->GetX();

  // reinitialize projectors
  iniProjectors();
  // reinitizize H0
  iniH0();

  // fetch material properties
  current_ye = body->ye(*track);
  current_density = body->density(*track);
  HI_constants = params.sqrt2*params.GF*params.Na*pow(params.cm,-3);

  // evolve projectors to current time
  EvolveProjectors(squids_time);

  // reading state
  H5LTget_dataset_info(group,"neustate", dims,nullptr,nullptr);
  std::unique_ptr<double[]> neudata(new double[dims[0]*dims[1]]);
  H5LTread_dataset_double(group,"neustate", neudata.get());

  H5LTget_dataset_info(group,"aneustate", dims,nullptr,nullptr);
  std::unique_ptr<double[]> aneudata(new double[dims[0]*dims[1]]);
  H5LTread_dataset_double(group,"aneustate", aneudata.get());

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
      size_t rounded_ne=round_up_to_aligned(ne);
      nc_factors.resize(std::vector<size_t>{nrhos,3,rounded_ne});
      if(tauregeneration){
        tau_hadlep_decays.resize(std::vector<size_t>{2,rounded_ne});
        tau_lep_decays.resize(std::vector<size_t>{2,rounded_ne});
      }
      if(iglashow)
        gr_factors.resize(std::vector<size_t>{rounded_ne});
      interactions_initialized = true;
    }
  }

  // read from user parameters
  hid_t user_parameters_id = H5Gopen(group, "user_parameters", H5P_DEFAULT);
  //H5Eset_auto (H5E_DEFAULT,(H5E_auto_t) H5Eprint,stderr);
  AddToReadHDF5(user_parameters_id);
  //H5Eset_auto (H5E_DEFAULT,nullptr,nullptr);
  H5Gclose(user_parameters_id);

  // we assume that this was created with the writer and got to this point!
  istate = true;
  ienergy = true;
  itrack = true;
  ibody = true;
  // initialize H0
  iniH0();
}

void nuSQUIDS::ReadStateHDF5(std::string str,std::string grp,std::string cross_section_grp_loc){
  //hid_t file_id,group_id,root_id;
  // open HDF5 file
  //std::cout << "reading from hdf5 file" << std::endl;
  H5File file(H5Fopen(str.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if(file < 0)
      throw std::runtime_error("nuSQUIDS::Error::file not found : " + str + ".");
  H5Handle rootGroup(H5Gopen(file, "/", H5P_DEFAULT), H5Gclose, "open HDF5 group");
  if (strncmp(grp.c_str(),"/",1)!=0){
    std::cout << "nuSQUIDS::ReadStateHDF5::Warning::group location name did not start with '/'. '/' was prepend to " + grp << std::endl;
    grp = "/"+grp;
  }
  H5Handle group(H5Gopen(rootGroup, grp.c_str(), H5P_DEFAULT),H5Gclose,"open HDF5 group");
  if( group < 0 )
      throw std::runtime_error("nuSQUIDS::Error::Group '" + grp + "' does not exist in HDF5.");

  // read number of neutrinos
  H5LTget_attribute_uint(group, "basic", "numneu", static_cast<unsigned int*>(&numneu));
  // neutrino/antineutrino/both
  int auxint;
  char auxchar[20];
  herr_t err;
  H5LTget_attribute_int(group, "basic", "NT", &auxint);
  NT = static_cast<NeutrinoType>(auxint);
  // interactions
  H5LTget_attribute_string(group,"basic","interactions", auxchar);
  std::string aux = auxchar;
  if ( aux == "True")
    iinteraction = true;
  else
    iinteraction = false;
  
  err=H5LTget_attribute_int(group, "basic", "oscillations", &auxint);
  if(err>=0) ioscillations=auxint;
  err=H5LTget_attribute_int(group, "basic", "tau_regeneration", &auxint);
  if(err>=0) tauregeneration=auxint;
  err=H5LTget_attribute_int(group, "basic", "glashow_resonance", &auxint);
  if(err>=0) iglashow=auxint;

  double squids_time;
  H5LTget_attribute_double(group, "basic", "squids_time", &squids_time);

  double squids_time_initial;
  H5LTget_attribute_double(group, "basic", "squids_time_initial", &squids_time_initial);

  // check version numbers
  unsigned int squids_version;
  H5LTget_attribute_uint(group, "basic", "squids_version_number", &squids_version);

  if ( squids_version > SQUIDS_VERSION )
    throw std::runtime_error("nuSQUIDS::ReadStateHDF5::Error: File was written using SQuIDS version " +
        std::to_string(squids_version) + " current version is " + std::to_string(SQUIDS_VERSION));

  unsigned int nusquids_version;
  H5LTget_attribute_uint(group, "basic", "nusquids_version_number", &nusquids_version);
  if ( nusquids_version > NUSQUIDS_VERSION )
    throw std::runtime_error("nuSQUIDS::ReadStateHDF5::Error: File was written using nuSQuIDS version " +
        std::to_string(nusquids_version) + " current version is " + std::to_string(NUSQUIDS_VERSION));

  // reading body and track
  if(nusquids_version>100000){
    H5Handle body_group(H5Gopen(group, "body", H5P_DEFAULT), H5Gclose, "open body group");
    std::string body_name = getStringAttribute(group,"body","name");
    body=GetBodyDeserializer(body_name)(body_group);

    H5Handle track_group(H5Gopen(group, "track", H5P_DEFAULT), H5Gclose, "open track group");
    std::string track_name = getStringAttribute(group,"track","name");
    track=GetTrackDeserializer(track_name)(track_group);
  } else {
    unsigned int body_id;
    hsize_t dimbody[1];
    H5LTget_attribute_uint(group,"body","ID",&body_id);

    H5LTget_dataset_info(group,"body",dimbody,nullptr,nullptr);
    double body_params[dimbody[0]];
    H5LTread_dataset_double(group,"body", body_params);

    hsize_t dimtrack[1];
    H5LTget_dataset_info(group,"track",dimtrack,nullptr,nullptr);
    std::unique_ptr<double[]> track_params(new double[dimtrack[0]]);
    H5LTread_dataset_double(group,"track", track_params.get());

    double x_current;
    H5LTget_attribute_double(group,"track","X",&x_current);

    // setting body and track
    SetBodyTrack(body_id,dimbody[0],body_params,dimtrack[0],track_params.get());

    // set trayectory to current time
    track->SetX(x_current);
  }

  // reading energy
  hsize_t dims[2];
  H5LTget_dataset_info(group, "energies", dims, nullptr, nullptr);

  ne = static_cast<unsigned int>(dims[0]);
  std::unique_ptr<double[]> energy_data(new double[ne]);
  H5LTread_dataset_double(group, "energies", energy_data.get());
  E_range = marray<double,1>{ne};
  for (unsigned int ie = 0; ie < ne; ie++)
    E_range[ie] = energy_data[ie];

  // initializing nuSQUIDS
  if (ne == 1){
    if(not inusquids)
      init(squids_time_initial);
    Set_E(energy_data[0]);
  }
  else {
    init(E_range,squids_time_initial);
  }

  // read and set mixing parameters
  for( unsigned int i = 0; i < numneu; i++ ){
    for( unsigned int j = i+1; j < numneu; j++ ){
      double th_value;
      std::string th_label = "th"+std::to_string(i+1)+std::to_string(j+1);
      H5LTget_attribute_double(group,"mixingangles", th_label.c_str(), &th_value);
      Set_MixingAngle(i,j,th_value);

      double delta_value;
      std::string delta_label = "delta"+std::to_string(i+1)+std::to_string(j+1);
      H5LTget_attribute_double(group,"CPphases", delta_label.c_str(), &delta_value);
      Set_CPPhase(i,j,delta_value);
    }
  }

  for( unsigned int i = 1; i < numneu; i++ ){
    double dm2_value;
    std::string dm2_label = "dm"+std::to_string(i+1)+"1sq";
    H5LTget_attribute_double(group,"massdifferences", dm2_label.c_str(), &dm2_value);
    Set_SquareMassDifference(i, dm2_value);
  }

  // reset current squids time
  Set_t(squids_time);
  // set time offset
  time_offset = squids_time - track->GetX();

  // reinitialize projectors
  iniProjectors();
  // reinitizize H0
  iniH0();

  // fetch material properties
  current_ye = body->ye(*track);
  current_density = body->density(*track);
  HI_constants = params.sqrt2*params.GF*params.Na*pow(params.cm,-3);

  // evolve projectors to current time
  EvolveProjectors(squids_time);

  // reading state
  H5LTget_dataset_info(group,"neustate", dims,nullptr,nullptr);
  std::unique_ptr<double[]> neudata(new double[dims[0]*dims[1]]);
  H5LTread_dataset_double(group,"neustate", neudata.get());

  H5LTget_dataset_info(group,"aneustate", dims,nullptr,nullptr);
  std::unique_ptr<double[]> aneudata(new double[dims[0]*dims[1]]);
  H5LTread_dataset_double(group,"aneustate", aneudata.get());

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
    H5Handle xs_grp;
    if ( cross_section_grp_loc == "") {
      xs_grp = H5Handle(H5Gopen(group, "crosssections", H5P_DEFAULT),H5Gclose,"open HDF5 crosssection group");
    } else {
      xs_grp = H5Handle(H5Gopen(rootGroup, cross_section_grp_loc.c_str(), H5P_DEFAULT),H5Gclose,"open HDF5 shared crosssection group");
    }
    
    int_struct = std::make_shared<InteractionStructure>();
	  
    // figure out how many target types we're dealing with
    unsigned int nTargets=0;
    std::vector<hsize_t> XSdim, dXSdim; //sizes to be determined
    if(nusquids_version<101100){
      nTargets=1; //old versions only understood isoscalar crosssections
      int_struct->targets={isoscalar_nucleon};
      
      XSdim.resize(3); //total xs tables were always rank 3
      dXSdim.resize(4); //differential xs tables were rank 4
    }
    else{
      //new files record the specific targets tabulated
      hsize_t rawNTargets=0;
      H5LTget_dataset_info(xs_grp,"targets",&rawNTargets,nullptr,nullptr);
      if(rawNTargets<1)
        throw std::runtime_error("Number of target types for cross sections does not make sense");
      nTargets=rawNTargets;
      int_struct->targets.resize(nTargets);
      H5Handle targetsDSet(H5Dopen2(xs_grp,"targets",H5P_DATASET_ACCESS_DEFAULT),H5Dclose,"open targets dataset");
      static_assert(sizeof(std::underlying_type<PDGCode>::type)==(32/CHAR_BIT), "Following code assumes PDGCode is int32_t");
      static_assert(std::is_signed<std::underlying_type<PDGCode>::type>::value, "Following code assumes PDGCode is int32_t");
      herr_t err=H5Dread(targetsDSet, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->targets.begin());
      if(err<0)
        throw std::runtime_error("Failed to read data from HDF5 dataset (targets)");
      
      XSdim.resize(4); //newer total xs tables are rank 4
      dXSdim.resize(5); //and differntial xs tables are rank 5
    }
    
    //get the actual dimensions of the main tables
    H5LTget_dataset_info(xs_grp,"sigmacc",&XSdim[0],nullptr,nullptr);
    H5LTget_dataset_info(xs_grp,"dNdEcc",&dXSdim[0],nullptr,nullptr);
    
    if(!nTargets) //if we don't yet know this, we're reading a new file and it's in the dataset dimensions
      nTargets=XSdim.front(); //first dimension is target type

    // initialize vectors
    InitializeInteractionVectors(nTargets);

    { // total crosssections
      H5Handle dsetSigmaCC(H5Dopen2(xs_grp,"sigmacc",H5P_DATASET_ACCESS_DEFAULT),H5Dclose,"open sigmacc dataset");
      H5Handle dsetSigmaNC(H5Dopen2(xs_grp,"sigmanc",H5P_DATASET_ACCESS_DEFAULT),H5Dclose,"open sigmanc dataset");
      H5Handle dsetSigmaGR(H5Dopen2(xs_grp,"sigmagr",H5P_DATASET_ACCESS_DEFAULT),H5Dclose,"open sigmagr dataset");
      
      //check the last dimension of the cross section arrays for padding
      if(int_struct->sigma_CC.extent(3)==ne){
        //no padding, which makes reading easy; dataspaces are the same shape
        //read sigma CC
        herr_t err=H5Dread(dsetSigmaCC, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->sigma_CC.begin());
        if(err<0)
          throw std::runtime_error("Failed to read data from HDF5 dataset (sigma CC)");
        
        //read sigma NC
        err=H5Dread(dsetSigmaNC, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->sigma_NC.begin());
        if(err<0)
          throw std::runtime_error("Failed to read data from HDF5 dataset (sigma NC)");
        
        //read sigma NC
        err=H5Dread(dsetSigmaGR, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->sigma_GR.begin());
        if(err<0)
          throw std::runtime_error("Failed to read data from HDF5 dataset (sigma GR)");
      }
      else{
        //After each run of data in the last dimension, the in-memory array has padding entries for 
        //alignment which do not exist on disk. Instead, we read into the blocks which fit between  
        //the sections of padding. 
        
        //The in-memory dataspace is a different size. 
        hsize_t mem_dim[4] {int_struct->sigma_CC.extent(0),
          int_struct->sigma_CC.extent(1),
          int_struct->sigma_CC.extent(2),
          int_struct->sigma_CC.extent(3)};
        H5Handle mem_dspace(H5Screate_simple(4,mem_dim,nullptr), H5Sclose, "create HDF5 dataspace");
        H5Handle mem_dspaceGR(H5Screate_simple(1,mem_dim+3,nullptr), H5Sclose, "create HDF5 dataspace");
        
        hsize_t stride[4]={1,1,1,1}; //we want everything, so all strides are minimal
        hsize_t count[4]={nTargets,nrhos,numneu,1}; //step through all but the last dimension taking blocks
        hsize_t block[4]={1,1,1,ne}; //the individual blocks are this size
        hsize_t start[4]={0,0,0,0};
        //designate which section of the memory we will write to, and the file contents we will write
        herr_t err=H5Sselect_hyperslab(mem_dspace, H5S_SELECT_SET, start, stride, count, block);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 in-memory dataspace");
        err=H5Sselect_hyperslab(mem_dspaceGR, H5S_SELECT_SET, start+3, stride+3, count+3, block+3);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 in-memory dataspace");
        
        err=H5Dread(dsetSigmaCC, H5T_NATIVE_DOUBLE, mem_dspace, H5S_ALL, H5P_DEFAULT, &*int_struct->sigma_CC.begin());
        if(err<0)
          throw std::runtime_error("Failed to read masked data from HDF5 dataset (sigma CC)");
        err=H5Dread(dsetSigmaNC, H5T_NATIVE_DOUBLE, mem_dspace, H5S_ALL, H5P_DEFAULT, &*int_struct->sigma_NC.begin());
        if(err<0)
          throw std::runtime_error("Failed to read masked data from HDF5 dataset (sigma NC)");
        err=H5Dread(dsetSigmaGR, H5T_NATIVE_DOUBLE, mem_dspaceGR, H5S_ALL, H5P_DEFAULT, &*int_struct->sigma_GR.begin());
        if(err<0)
          throw std::runtime_error("Failed to read masked data from HDF5 dataset (sigma GR)");
      }
    }

    { // differential crosssections
      H5Handle dsetDSigmaCC(H5Dopen2(xs_grp,"dNdEcc",H5P_DATASET_ACCESS_DEFAULT),H5Dclose,"open dNdEcc dataset");
      H5Handle dsetDSigmaNC(H5Dopen2(xs_grp,"dNdEnc",H5P_DATASET_ACCESS_DEFAULT),H5Dclose,"open dNdEnc dataset");
      H5Handle dsetDSigmaGR(H5Dopen2(xs_grp,"dNdEgr",H5P_DATASET_ACCESS_DEFAULT),H5Dclose,"open dNdEgr dataset");
      
      //check for padding
      if(int_struct->dNdE_CC.extent(4)==ne){
        //no padding
        herr_t err=H5Dread(dsetDSigmaCC, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_CC.begin());
        if(err<0)
          throw std::runtime_error("Failed to read data from HDF5 dataset (dNdE CC)");
        
        err=H5Dread(dsetDSigmaNC, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_NC.begin());
        if(err<0)
          throw std::runtime_error("Failed to read data from HDF5 dataset (dNdE NC)");
        
        err=H5Dread(dsetDSigmaGR, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_GR.begin());
        if(err<0)
          throw std::runtime_error("Failed to read data from HDF5 dataset (dNdE GR)");
      }
      else{ //need to avoid padding
        //in-memory dataspaces 
        hsize_t mem_dim[5] {int_struct->dNdE_CC.extent(0),
          int_struct->dNdE_CC.extent(1),
          int_struct->dNdE_CC.extent(2),
          int_struct->dNdE_CC.extent(3),
          int_struct->dNdE_CC.extent(4)};
        H5Handle mem_dspace(H5Screate_simple(5,mem_dim,nullptr), H5Sclose, "create HDF5 dataspace");
        H5Handle mem_dspaceGR(H5Screate_simple(2,mem_dim+3,nullptr), H5Sclose, "create HDF5 dataspace");
        
        hsize_t stride[5]={1,1,1,1,1}; //we want everything, so all strides are minimal
        hsize_t count[5]={nTargets,nrhos,numneu,ne,1}; //step through all but the last dimension taking blocks
        hsize_t block[5]={1,1,1,1,ne}; //the individual blocks are this size
        hsize_t start[5]={0,0,0,0,0};
        //designate which section of the memory we will write to, and the file contents we will write
        herr_t err=H5Sselect_hyperslab(mem_dspace, H5S_SELECT_SET, start, stride, count, block);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 in-memory dataspace");
        err=H5Sselect_hyperslab(mem_dspaceGR, H5S_SELECT_SET, start+3, stride+3, count+3, block+3);
        if(err<0)
          throw std::runtime_error("Failed to set selection in HDF5 in-memory dataspace");
        
        err=H5Dread(dsetDSigmaCC, H5T_NATIVE_DOUBLE, mem_dspace, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_CC.begin());
        if(err<0)
          throw std::runtime_error("Failed to read masked data from HDF5 dataset (sigma CC)");
        err=H5Dread(dsetDSigmaNC, H5T_NATIVE_DOUBLE, mem_dspace, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_NC.begin());
        if(err<0)
          throw std::runtime_error("Failed to read masked data from HDF5 dataset (sigma NC)");
        err=H5Dread(dsetDSigmaGR, H5T_NATIVE_DOUBLE, mem_dspaceGR, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_GR.begin());
        if(err<0)
          throw std::runtime_error("Failed to read masked data from HDF5 dataset (sigma GR)");
      }
    }
    
    { // tau decay data
      H5Handle dNdEtauall(H5Dopen2(xs_grp,"dNdEtauall",H5P_DATASET_ACCESS_DEFAULT),H5Dclose,"open dNdEtauall dataset");
      H5Handle dNdEtaulep(H5Dopen2(xs_grp,"dNdEtaulep",H5P_DATASET_ACCESS_DEFAULT),H5Dclose,"open dNdEtaulep dataset");
      
      int tau_spectra_rank;
      H5LTget_dataset_ndims(xs_grp,"dNdEtauall",&tau_spectra_rank);
      
      if(tau_spectra_rank == 3){ // spectra with polarization
        //check for padding
        if(int_struct->dNdE_tau_all.extent(2)==ne){
          //no padding
          herr_t err=H5Dread(dNdEtauall, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_tau_all.begin());
          if(err<0)
            throw std::runtime_error("Failed to read data from HDF5 dataset (dNdE tau all)");
          
          err=H5Dread(dNdEtaulep, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_tau_lep.begin());
          if(err<0)
            throw std::runtime_error("Failed to read data from HDF5 dataset (dNdE tau leptonic)");
        }
        else{ //harder case with padding
          //in-memory dataspace 
          hsize_t source_dim[3] {int_struct->dNdE_tau_all.extent(0),
            int_struct->dNdE_tau_all.extent(1),
            int_struct->dNdE_tau_all.extent(2)};
          H5Handle mem_dspace(H5Screate_simple(3,source_dim,nullptr), H5Sclose, "create HDF5 dataspace");
          
          hsize_t stride[3]={1,1,1}; //we want everything, so all strides are minimal
          hsize_t count[3]={nrhos,ne,1}; //step through all but the last dimension taking blocks
          hsize_t block[3]={1,1,ne}; //the individual blocks are this size
          hsize_t start[3]={0,0,0};
          //designate which section of the memory we will read, and the file contents we will write
          herr_t err=H5Sselect_hyperslab(mem_dspace, H5S_SELECT_SET, start, stride, count, block);
          if(err<0)
            throw std::runtime_error("Failed to set selection in HDF5 in-memory dataspace");
          
          err=H5Dread(dNdEtauall, H5T_NATIVE_DOUBLE, mem_dspace, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_tau_all.begin());
          if(err<0)
            throw std::runtime_error("Failed to read masked data from HDF5 dataset (dNdE tau all)");
          err=H5Dread(dNdEtaulep, H5T_NATIVE_DOUBLE, mem_dspace, H5S_ALL, H5P_DEFAULT, &*int_struct->dNdE_tau_lep.begin());
          if(err<0)
            throw std::runtime_error("Failed to read masked data from HDF5 dataset (dNdE tau leptonic)");
        }
      }
      else if (tau_spectra_rank == 2) { // legacy data with polarization averaged
        hsize_t dNdEtaudim[2];
        H5LTget_dataset_info(xs_grp,"dNdEtauall", dNdEtaudim,nullptr,nullptr);
        
        std::unique_ptr<double[]> dNdEtauall(new double[dNdEtaudim[0]*dNdEtaudim[1]]);
        H5LTread_dataset_double(xs_grp,"dNdEtauall", dNdEtauall.get());
        std::unique_ptr<double[]> dNdEtaulep(new double[dNdEtaudim[0]*dNdEtaudim[1]]);
        H5LTread_dataset_double(xs_grp,"dNdEtaulep", dNdEtaulep.get());
        
        for(unsigned int rho = 0; rho < nrhos; rho++){
          for( unsigned int e1 = 0; e1 < ne; e1++){
            for( unsigned int e2 = 0; e2 < e1; e2++){
              int_struct->dNdE_tau_all[rho][e1][e2] = dNdEtauall[e1*ne + e2];
              int_struct->dNdE_tau_lep[rho][e1][e2] = dNdEtaulep[e1*ne + e2];
            }
          }
        }
      }
      else {
        throw std::runtime_error("nuSQUIDS::ReadStateHDF5: Error. Malformed HDF5 files when trying to deserialize tau decay distributions.");
      }
    }

    nc_factors.resize(std::vector<size_t>{nrhos,3,ne});
    if(tauregeneration){
      tau_hadlep_decays.resize(std::vector<size_t>{2,ne});
      tau_lep_decays.resize(std::vector<size_t>{2,ne});
    }
    if(iglashow)
      gr_factors.resize(std::vector<size_t>{ne});

    interactions_initialized = true;
  }

  // read from user parameters
  hid_t user_parameters_id = H5Gopen(group, "user_parameters", H5P_DEFAULT);
  //H5Eset_auto (H5E_DEFAULT,(H5E_auto_t) H5Eprint,stderr);
  AddToReadHDF5(user_parameters_id);
  //H5Eset_auto (H5E_DEFAULT,nullptr,nullptr);
  H5Gclose(user_parameters_id);

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
          auto earth = std::make_shared<EarthAtm>();
          body = earth;
          // track_param[2] corresponds to the zenith angle
          track = std::make_shared<EarthAtm::Track>(earth->MakeTrack(track_params[2]));
          break;
        }
      default:
        {
          std::cerr << "nuSQUIDS::SetBodyTrack : unknown body/track" << std::endl;
          exit(1);
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
  if( not iinteraction and opt )
    throw std::runtime_error("nuSQUIDS::Error::nuSQuIDs has been initialized without interactions, thus tau regeneration cannot be enabled.");
  if( numneu < 3 and opt )
    throw std::runtime_error("nuSQUIDS::Error::nuSQuIDs at least three neutrino flavors must be included to enable tau regeneration.");
  tauregeneration = opt;
}
  
void nuSQUIDS::Set_GlashowResonance(bool opt){
  if( not iinteraction and opt )
    throw std::runtime_error("nuSQUIDS::Error::nuSQuIDs has been initialized without interactions, thus the Glashow resonance cannot be enabled.");
  iglashow = opt;
}

void nuSQUIDS::Set_ProgressBar(bool opt){
  progressbar = opt;
}

void nuSQUIDS::Set_EvolLowPassCutoff(double val){
  evol_lowpass_cutoff = val;
}

void nuSQUIDS::Set_EvolLowPassScale(double val){
  evol_lowpass_scale = val;
}
void nuSQUIDS::Set_IncludeOscillations(bool opt){
  ioscillations = opt;
}

void nuSQUIDS::Set_AllowConstantDensityOscillationOnlyEvolution(bool opt){
  allowConstantDensityOscillationOnlyEvolution = opt;
}

std::shared_ptr<Track> nuSQUIDS::GetTrack(){
  return track;
}

std::shared_ptr<Body> nuSQUIDS::GetBody(){
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
    throw std::invalid_argument("nuSQUIDS::Set_SquareMassDifference::Error: Index greater than number of neutrino flavors.");
  params.SetEnergyDifference(i,val);
  istate = false;
}

double nuSQUIDS::Get_SquareMassDifference( unsigned int i ) const {
  if ( i > numneu )
    throw std::invalid_argument("nuSQUIDS::Set_SquareMassDifference::Error: Index greater than number of neutrino flavors.");
  return params.GetEnergyDifference(i);
}

void nuSQUIDS::Set_MixingParametersToDefault(void){
  // set parameters as in NuFit 5.0 (arXiv:2007.14792) NO
  // but with delta_CP = 0.0
  Set_MixingAngle(0,1,0.583638); // th12
  Set_MixingAngle(0,2,0.149575); // th13
  Set_MixingAngle(1,2,0.855211); // th23

  Set_SquareMassDifference(1,7.420e-05); // dm^2_21
  Set_SquareMassDifference(2,0.002514); // dm^2_31

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
NT(other.NT),
int_struct(std::move(other.int_struct)),
int_state(std::move(other.int_state)),
positivization_scale(other.positivization_scale),
body(other.body),
track(other.track),
DM2(other.DM2),
H0_array(std::move(other.H0_array)),
HI_constants(other.HI_constants),
current_density(other.current_density),
current_ye(other.current_ye),
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
ioscillations(other.ioscillations),
tauregeneration(other.tauregeneration),
iglashow(other.iglashow),
positivization(other.positivization),
progressbar(other.progressbar),
progressbar_count(other.progressbar_count),
progressbar_loop(other.progressbar_loop)
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
  int_state = std::move(other.int_state);
  positivization_scale = other.positivization_scale;
  body = other.body;
  track = other.track;
  DM2 = other.DM2;
  H0_array = std::move(other.H0_array);
  HI_constants = other.HI_constants;
  current_density = other.current_density;
  current_ye = other.current_ye;
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
  ioscillations = other.ioscillations;
  tauregeneration = other.tauregeneration;
  iglashow = other.iglashow;
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
