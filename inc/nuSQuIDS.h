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


#ifndef __nuSQUID_H
#define __nuSQUID_H

#if __cplusplus < 201103L
#error C++11 compiler required. Update your compiler and use the flag -std=c++11
#endif

#include "version.h"
#include "body.h"
#include "xsections.h"
#include "taudecay.h"
#include "marray.h"
#include "aligned_alloc.h"

#include <algorithm>
#include <cstring>
#include <SQuIDS/SQuIDS.h>
#include <memory>
#include <map>
#include <stdexcept>

#include "H5Epublic.h"
#include "H5Tpublic.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Gpublic.h"
#include "H5Fpublic.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define ManualTauReinjection
//#define FixCrossSections

//#define UpdateInteractions_DEBUG

namespace nusquids{

enum NeutrinoType {
  neutrino=0b01,
  antineutrino=0b10,
  both=0b11
};

enum Basis {
  mass=0b01,
  flavor=0b10,
  interaction=0b11
};


///\class nuSQUIDS
///\brief nu-SQuIDS main class
class nuSQUIDS: public squids::SQuIDS {
  // nuSQUIDSAtm is a friend class so it can use H0 to evolve operators
  // and thus evaluate expectation values.
  template<typename,typename>
  friend class nuSQUIDSAtm;
protected:
  /// \brief Sets the basis in which the problem will be solved.
  ///
  /// If interaction basis is used the projectors will be evolved at
  /// every time step. On the other hand, if mass basis is used no evolution
  /// is performed.
    Basis basis = interaction;
    /// \brief number of neutrino flavors.
    unsigned int numneu;
    /// \brief number of energy nodes.
    unsigned int ne = nx;

    /// \brief Returns the number of nucleons at a given position.
    ///
    /// Isoscalar medium is assumed and the position is obtained
    /// at x = track.GetX().
    double GetNucleonNumber() const;

    /// \brief Updates the interaction length arrays.
    ///
    /// Uses GetNucleonNumber() together with the stored cross section
    /// information to update: nuSQUIDS#invlen_NC, nuSQUIDS#invlen_CC, nuSQUIDS#invlen_INT, and nuSQUIDS#invlen_tau.
    void UpdateInteractions();

    /// \brief Contains the energy nodes.
    marray<double,1> E_range;

    /// \brief Contains the spaces between nodes.
    ///
    /// It has length len(E_range)-1.
    marray<double,1> delE;

    /// \brief Interface that calculate and interpolates neutrino cross sections.
    std::shared_ptr<NeutrinoCrossSections> ncs;

    /// \brief Interface that calculate and interpolates tau decay spectral functions.
    TauDecaySpectra tdc;

    /// \brief NT keeps track if the problem consists of neutrinos, antineutrinos, or both.
    NeutrinoType NT = both;
  
    ///aligned arrays should be aligned to at least this many items
    constexpr static size_t preferred_alignment=4;
    constexpr static size_t log2(size_t n) {
      return((n>1) ? log2(n/2)+1 : 0);
    }
    constexpr static size_t round_up_to_aligned(size_t n){
      return((n&~(preferred_alignment-1))+preferred_alignment);
    }

  public:
    /// \brief Struct that contains all cross section information.
    struct InteractionStructure {
        /// \brief Neutrino charge current differential cross section with respect to
        /// the outgoing lepton energy.
        ///
        /// The array four dimensional is constructed when InitializeInteractionVectors() is called and
        /// its initialized when InitializeInteractions() is called. The first dimension
        /// is number of neutrino types (neutrino/antineutrino/both), the second the neutrino flavor,
        /// and the last two the initial and final energy node respectively.
        marray<double,4,aligned_allocator<double>> dNdE_CC;
        /// \brief Neutrino neutral current differential cross section with respect to
        /// the outgoing lepton energy.
        ///
        /// The array four dimensional is constructed when InitializeInteractionVectors() is called and
        /// its initialized when InitializeInteractions() is called. The first dimension
        /// is number of neutrino types (neutrino/antineutrino/both), the second the neutrino flavor,
        /// and the last two the initial and final energy node respectively.
        marray<double,4,aligned_allocator<double>> dNdE_NC;
        /// \brief Glashow resonance differential cross section (electron antineutrino only
        /// \details Dimensions are initial and final energy node
        marray<double,2,aligned_allocator<double>> dNdE_GR;
        /// \brief Array that contains the inverse of the neutrino charge current mean free path.
        /// \details The array contents are in natural units (i.e. eV) and is update when
        /// UpdateInteractions() is called. The first dimension corresponds to the neutrino type,
        /// the second to the flavor, and the last one to the energy.
        marray<double,3,aligned_allocator<double>> invlen_CC;
        /// \brief Array that contains the inverse of the neutrino neutral current mean free path.
        /// \details The array contents are in natural units (i.e. eV) and is update when
        /// UpdateInteractions() is called. The first dimension corresponds to the neutrino type,
        /// the second to the flavor, and the last one to the energy.
        marray<double,3,aligned_allocator<double>> invlen_NC;
        /// \brief Glashow resonance inverse interaction length (electron antineutrino only)
        /// \details 1 entry per energy node
        marray<double,1,aligned_allocator<double>> invlen_GR;
        /// \brief Array that contains the inverse of the neutrino total mean free path.
        /// \details The array contents are in natural units (i.e. eV) and is update when
        /// UpdateInteractions() is called. Numerically it is just nuSQUIDS::invlen_NC and nuSQUIDS::invlen_CC
        /// added together.
        marray<double,3,aligned_allocator<double>> invlen_INT;
        /// \brief Array that contains the neutrino charge current cross section.
        /// \details The first dimension corresponds to the neutrino type, the second to the flavor, and
        /// the final one to the energy node. Its contents are in natural units, i.e. eV^-2. It is
        /// initialized by InitializeInteractions() .
        marray<double,3,aligned_allocator<double>> sigma_CC;
        /// \brief Array that contains the neutrino neutral current cross section.
        /// \details The first dimension corresponds to the neutrino type, the second to the flavor, and
        /// the final one to the energy node. Its contents are in natural units, i.e. eV^-2. It is
        /// initialized by InitializeInteractions() .
        marray<double,3,aligned_allocator<double>> sigma_NC;
        /// \brief Glashow resonance cross section (electron antineutrino only)
        /// \details 1 entry per energy node, in natural units
        marray<double,1,aligned_allocator<double>> sigma_GR;
        /// \brief Array that contains the inverse of the tau decay length for each energy node.
        marray<double,1,aligned_allocator<double>> invlen_tau;
        /// \brief Array that contains the tau decay spectrum to all particles.
        /// \details The first dimension corresponds to initial tau energy and the
        /// second one to the outgoing lepton.
        marray<double,2,aligned_allocator<double>> dNdE_tau_all;
        /// \brief Array that contains the tau decay spectrum to leptons.
        /// \brief Array that contains the tau decay spectrum to all particles.
        /// \details The first dimension corresponds to initial tau energy and the
        /// second one to the outgoing lepton.
        marray<double,2,aligned_allocator<double>> dNdE_tau_lep;
        /// \brief Default constructor
        InteractionStructure():
        dNdE_CC(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        dNdE_NC(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        dNdE_GR(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        invlen_CC(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        invlen_NC(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        invlen_GR(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        invlen_INT(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        sigma_CC(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        sigma_NC(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        sigma_GR(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        invlen_tau(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        dNdE_tau_all(aligned_allocator<double>(log2(preferred_alignment*sizeof(double)))),
        dNdE_tau_lep(aligned_allocator<double>(log2(preferred_alignment*sizeof(double))))
        {}
        /// \brief Assignment operator
        InteractionStructure& operator=(const InteractionStructure& other){
          dNdE_CC=other.dNdE_CC;
          dNdE_NC=other.dNdE_NC;
          dNdE_GR=other.dNdE_GR;
          invlen_NC=other.invlen_NC;
          invlen_CC=other.invlen_CC;
          invlen_GR=other.invlen_GR;
          invlen_INT=other.invlen_INT;
          sigma_CC=other.sigma_CC;
          sigma_NC=other.sigma_NC;
          sigma_GR=other.sigma_GR;
          invlen_tau=other.invlen_tau;
          dNdE_tau_all=other.dNdE_tau_all;
          dNdE_tau_lep=other.dNdE_tau_lep;
          return(*this);
        }
        /// \brief Move assignment operator
        InteractionStructure& operator=(InteractionStructure&& other){
          dNdE_CC=std::move(other.dNdE_CC);
          dNdE_NC=std::move(other.dNdE_NC);
          dNdE_GR=std::move(other.dNdE_GR);
          invlen_NC=std::move(other.invlen_NC);
          invlen_CC=std::move(other.invlen_CC);
          invlen_GR=std::move(other.invlen_GR);
          invlen_INT=std::move(other.invlen_INT);
          sigma_CC=std::move(other.sigma_CC);
          sigma_NC=std::move(other.sigma_NC);
          sigma_GR=std::move(other.sigma_GR);
          invlen_tau=std::move(other.invlen_tau);
          dNdE_tau_all=std::move(other.dNdE_tau_all);
          dNdE_tau_lep=std::move(other.dNdE_tau_lep);
          return(*this);
        }
        /// \brief Equality operator
        bool operator==(const InteractionStructure& other) const{
          if (
               dNdE_CC.size() != other.dNdE_CC.size() or
               dNdE_NC.size() != other.dNdE_NC.size() or
               dNdE_GR.size() != other.dNdE_GR.size() or
               invlen_CC.size() != other.invlen_CC.size() or
               invlen_NC.size() != other.invlen_NC.size() or
               invlen_GR.size() != other.invlen_GR.size() or
               invlen_INT.size() != other.invlen_INT.size() or
               sigma_CC.size() != other.sigma_CC.size() or
               sigma_NC.size() != other.sigma_NC.size() or
               sigma_GR.size() != other.sigma_GR.size() or
               invlen_tau.size() != other.invlen_tau.size() or
               dNdE_tau_all.size() != other.dNdE_tau_all.size() or
               dNdE_tau_lep.size() != other.dNdE_tau_lep.size()
             )
          {
            return false;
          }

          if ( not std::equal(dNdE_CC.begin(),dNdE_CC.end(),other.dNdE_CC.begin()) )
            return false;

          if ( not std::equal(dNdE_NC.begin(),dNdE_NC.end(),other.dNdE_NC.begin()) )
            return false;

          if ( not std::equal(dNdE_GR.begin(),dNdE_GR.end(),other.dNdE_GR.begin()) )
            return false;

          if ( not std::equal(invlen_CC.begin(),invlen_CC.end(),other.invlen_CC.begin()) )
            return false;

          if ( not std::equal(invlen_NC.begin(),invlen_NC.end(),other.invlen_NC.begin()) )
            return false;

          if ( not std::equal(invlen_GR.begin(),invlen_GR.end(),other.invlen_GR.begin()) )
            return false;

          if ( not std::equal(invlen_INT.begin(),invlen_INT.end(),other.invlen_INT.begin()) )
            return false;

          if ( not std::equal(sigma_CC.begin(),sigma_CC.end(),other.sigma_CC.begin()) )
            return false;

          if ( not std::equal(sigma_NC.begin(),sigma_NC.end(),other.sigma_NC.begin()) )
            return false;

          if ( not std::equal(sigma_GR.begin(),sigma_GR.end(),other.sigma_GR.begin()) )
            return false;

          if ( not std::equal(invlen_tau.begin(),invlen_tau.end(),other.invlen_tau.begin()) )
            return false;

          if ( not std::equal(dNdE_tau_all.begin(),dNdE_tau_all.end(),other.dNdE_tau_all.begin()) )
            return false;

          if ( not std::equal(dNdE_tau_lep.begin(),dNdE_tau_lep.end(),other.dNdE_tau_lep.begin()) )
            return false;
          // all is good - lets roll
          return true;
        }
        /// \brief Inequality operator
        bool operator!=(const InteractionStructure& other) const {
          return not (*this==other);
        }
        /// \brief Copy constructor operator
        InteractionStructure(InteractionStructure& other):
          dNdE_CC(other.dNdE_CC),dNdE_NC(other.dNdE_NC),dNdE_GR(other.dNdE_GR),
          invlen_CC(other.invlen_CC),invlen_NC(other.invlen_NC),invlen_GR(other.invlen_GR),invlen_INT(other.invlen_INT),
          sigma_CC(other.sigma_CC),sigma_NC(other.sigma_NC),sigma_GR(other.sigma_GR),
          invlen_tau(other.invlen_tau),
          dNdE_tau_all(other.dNdE_tau_all),
          dNdE_tau_lep(other.dNdE_tau_lep)
        {}
        /// \brief Move constructor operator
        InteractionStructure(InteractionStructure&& other):
          dNdE_CC(std::move(other.dNdE_CC)),dNdE_NC(std::move(other.dNdE_NC)),dNdE_GR(std::move(other.dNdE_GR)),
          invlen_CC(std::move(other.invlen_CC)),invlen_NC(std::move(other.invlen_NC)),invlen_GR(std::move(other.invlen_GR)),invlen_INT(std::move(other.invlen_INT)),
          sigma_CC(std::move(other.sigma_CC)),sigma_NC(std::move(other.sigma_NC)),sigma_GR(std::move(other.sigma_GR)),
          invlen_tau(std::move(other.invlen_tau)),
          dNdE_tau_all(std::move(other.dNdE_tau_all)),
          dNdE_tau_lep(std::move(other.dNdE_tau_lep))
        {}
    };
  protected:
    /// \brief Interaction struct that contains all cross section information tabulated.
    std::shared_ptr<InteractionStructure> int_struct;

    /// \brief Length upon which the neutrino fluxes will be positivized.
    double positivization_scale;

    /// \brief Body where the neutrino propagation takes place.
    std::shared_ptr<Body> body;
    /// \brief Trayectory within the body.
    /// \details Stores the position within the body and its updated every evolution
    /// step.
    std::shared_ptr<Track> track;

    /// \brief SU_vector that represents the neutrino square mass difference matrix in the mass basis.
    ///  It is used to construct nuSQUIDS#H0_array and H0()
    squids::SU_vector DM2;
    /// \brief Stores the time independent hamiltonian corresponding to each energy node.
    marray<squids::SU_vector,1> H0_array;
  
    /// \brief Product of physical constants needed by HI.
    double HI_constants;
    /// \brief The density of the body at the current point on the track
    double current_density;
    /// \brief The electron fraction of the body at the current point on the track
    double current_ye;

    /// \brief Mass basis projectors.
    /// \details The i-entry corresponds to the projector in the ith mass eigenstate.
    marray<squids::SU_vector,1> b0_proj;
    /// \brief Flavor basis projectors.
    /// \details The first dismension corresponds to the neutrino type. When NeutrinoType = both, then
    /// the first dimension equal 0 means neutrinos and 1 antineutrinos. The second dimension corresponds
    /// to the flavor eigenstates where 0 corresponds to nu_e, 1 to nu_mu, 2 to nu_tau, and the others
    /// to sterile neutrinos.
    marray<squids::SU_vector,2> b1_proj;
    /// \brief Mass basis projectors in the interaction picture.
    /// The index meaning are the same as nuSQUIDS#b1_proj but for mass eigenstates,
    /// with an added third dimension that corresponds to the energy node index.
    marray<squids::SU_vector,3> evol_b0_proj;
    /// \brief Flavor basis projectors in the interaction picture.
    /// The index meaning are the same as nuSQUIDS#b1_proj ,
    /// with an added third dimension that corresponds to the energy node index.
    marray<squids::SU_vector,3> evol_b1_proj;

    ///Cache of precalculated results for InteractionsRho.
    ///Used only when interactions are on and oscillations are off.
    marray<squids::SU_vector,2> interaction_cache;
    ///Backing storage for the SU_vectors in interaction_cache so that each can
    ///be ensured to be correctly aligned for best performance.
    std::unique_ptr<double[]> interaction_cache_store;
    ///Total size of interaction_cache_store
    size_t interaction_cache_store_size;
  
    ///Initialize the interaction cahce data structures to the corretc sizes for
    ///the current problem.
    ///Should be invoked only when interactions are on and oscillations are off.
    void SetUpInteractionCache();

    /// \brief Evolves the flavor projectors in the interaction basis to a time t.
    /// \details It uses H0() to evolve SQUIDS#b0_proj and SQUIDS#b1_proj into 
    /// SQUIDS#evol_b0_proj and SQUIDS#evol_b1_proj.
    /// \warning Since the RHS of the differential equation only involves flavor projectors
    /// we do not current evolve mass projectors.
    void EvolveProjectors(double t);

    // bool requirements
  private:
    /// \brief Boolean that signals that the full hamiltonian including extra potentials should be used for 
    /// the projector evolution
    bool use_full_hamiltonian_for_projector_evolution = false;
    /// \brief Boolean that toggles debug information to be printed
    bool debug = false;
    /// \brief Boolean that signals the object correct initialization.
    bool inusquids = false;
    /// \brief Boolean that signals that a Body object has being set.
    bool ibody = false;
    /// \brief Boolean that signals that the neutrino energies has being set.
    bool ienergy = false;
    /// \brief Boolean that signals that a Track object has being set.
    bool itrack = false;
    /// \brief Boolean that signals that an initial state has being set.
    bool istate = false;
    /// \brief Boolean that signals that interactions will be taken into account.
    bool iinteraction = false;
    /// \brief Whether interactions have been initialized
    bool interactions_initialized = false;
    /// \brief Boolean that signals that neutrino oscillations will be taken into account.
    bool ioscillations = true;
    /// \brief Boolean that signals that tau regeneration is being used.
    bool tauregeneration = false;
    /// \brief Whether the Glashsow resonance for electron antineutrinos is used.
    bool iglashow = false;
    /// \brief Boolean that signals that positivization will be enforced.
    bool positivization = false;
    /// \brief Boolean that signals that a progress bar will be printed.
    bool progressbar = false;
    /// \brief Integer to keep track of the progress bar evolution.
    int progressbar_count = 0;
    /// \brief Number of steps upon which the progress bar will be updated.
    int progressbar_loop = 100;
    /// \brief Time offset between SQuIDS time and Track(x).
    double time_offset;
    /// \brief Force flavor projections to be positive.
    void PositivizeFlavors();
  protected:
    /// \brief Initializes flavor and mass projectors
    /// \warning Antineutrinos are handle by means of the AntineutrinoCPFix() function
    /// which temporary flips the CP phase.
    void iniProjectors();
    /// \brief Reinitializes the flavor projectors.
    /// \details It is called before evolving the system and every time the system
    /// state is going to be set, since the definition of mass and flavor basis might
    /// have changed.
    void SetIniFlavorProyectors();
    /// \brief Initializes the time independent hamiltonian for each energy node.
    /// \details It constructs nuSQUIDS#DM2 and nuSQUIDS#H0_array
    void iniH0();
    /// \brief Changes the CP sign of the complex matrices stored in params.
    void AntineutrinoCPFix(unsigned int irho);
    /// \brief Serializes the initialization of the Body and Track objects.
    /// @see ReadStateHDF5
    void SetBodyTrack(unsigned int,unsigned int,double*,unsigned int,double*);

    /// \brief General initilizer for the multi energy mode
    /// @param E_vector Energy nodes [eV].
    /// @param xini The initial position of the system.
    void init(marray<double,1> E_vector,double xini=0.0);

    /// \brief Initilizer for the single energy mode
    /// @param xini The initial position of the system. By default is set to 0.
    void init(double xini=0.0);

    /// \brief Initilizes auxiliary cross section arrays.
    /// \details The arrays are initialize, but not filled with contented. To fill the arrays
    /// call InitializeInteractions().
    void InitializeInteractionVectors();

    /// \brief Fills in auxiliary cross section arrays.
    /// \details It uses nuSQUIDS#ncs and nuSQUIDS#tdc to fill in the values of
    /// nuSQUIDS#dNdE_CC , nuSQUIDS#dNdE_NC , nuSQUIDS#sigma_CC , nuSQUIDS#sigma_NC ,
    /// nuSQUIDS#invlen_CC , nuSQUIDS#invlen_NC , nuSQUIDS#invlen_INT ,
    /// nuSQUIDS#dNdE_tau_all , nuSQUIDS#dNdE_tau_lep , nuSQUIDS#invlen_tau
    /// in int_struct.
    /// @see InitializeInteractionVectors
    void GetCrossSections();
  private:

    /// \brief Prints progress bar.
    /// \details To enable it call Set_ProgressBar()
    void ProgressBar() const;
  protected:
    /// \brief Updates all quantities needed during the evaluation of the RHS.
    /// @param x Position of the system.
    /// \details Dependind on the basis used it evolves the flavor projectors
    /// by means of EvolveProjectors(). Also, when interactions are considered,
    /// it updates the interaction arrays by calling UpdateInteractions(). Finally,
    /// it handles control to the user implemented updates via AddToPreDerive().
    /// @see AddToPreDerive
    /// @see EvolveProjectors
    /// @see UpdateInteractions
    void PreDerive(double x);
    /// \brief User supplied function that is called before performing a derivative.
    /// @param x Position of the system.
    /// @see PreDerive
    virtual void AddToPreDerive(double x){}

    /// \brief Multiple energy mode constructor using precalculated cross section information.
    /// @param E_vector Energy nodes [eV].
    /// @param numneu Number of neutrino flavors.
    /// @param NT NeutrinoType: neutrino,antineutrino, or both (simultaneous solution).
    /// @param iinteraction Sets the neutrino noncoherent neutrino interactions on.
    /// @param int_struct Structure that contains cross section information.
    /// @warning The user should assure that the cross section struct has the same range as the provided energy range
    /// @see init
    /// \todo put asserts here to ensure some minimal safety 
    nuSQUIDS(marray<double,1> E_vector,unsigned int numneu,NeutrinoType NT,std::shared_ptr<InteractionStructure> int_struct,
             bool iinteraction = true):
    numneu(numneu),iinteraction(iinteraction),NT(NT),int_struct(int_struct)
    {
     // assert(int_struct.
      init(E_vector);
    }

    /// \brief Reads and constructs the object from an HDF5 file.
    /// @param hdf5_filename Filename of the HDF5 to use for construction.
    /// @param group Path to the group where the nuSQUIDS content will be saved.
    /// @param iis InteractionStructure that contains valid cross sections.
    /// \details By default all contents are assumed to be in the \c root of the HDF5 file
    /// if the user wants a different location it can use \c group and \c save_cross_sections to
    /// change the nuSQUIDS location and the cross section information location.
    /// Finally, it calls AddToReadHDF5() to enable user defined parameters
    /// to be read.
    /// @see AddToReadHDF5
    /// @see WriteStateHDF5
    void ReadStateHDF5(std::string hdf5_filename,std::string group = "/", std::shared_ptr<InteractionStructure> iis = nullptr);

  public:
    /************************************************************************************
     * CONSTRUCTORS
    *************************************************************************************/

    /// \brief Default void constructor.
    nuSQUIDS(){}

    /// \brief Move constructor.
    nuSQUIDS(nuSQUIDS&&);

    /// \brief Multiple energy mode constructor.
    /// @param E_vector Energy nodes [eV].
    /// @param numneu Number of neutrino flavors.
    /// @param NT NeutrinoType: neutrino,antineutrino, or both (simultaneous solution).
    /// @param iinteraction Sets the neutrino noncoherent neutrino interactions on.
    /// @param ncs Cross section object.
    /// \details By default the energy scale is logarithmic and interactions are turn off.
    /// \warning When interactions are present interpolation is performed to precalculate the neutrino
    /// cross section which make take considertable time depending on the energy grid.
    /// @see init
    nuSQUIDS(marray<double,1> E_vector,unsigned int numneu,NeutrinoType NT = both,
       bool iinteraction = false, std::shared_ptr<NeutrinoCrossSections> ncs = nullptr):
    numneu(numneu),ncs(ncs),iinteraction(iinteraction),NT(NT)
    {init(E_vector);}

    /// \brief Single energy mode constructor.
    /// @param numneu Number of neutrino flavors.
    /// @param NT NeutrinoType: neutrino or antineutrino.
    /// \warning Interactions are not possible in the single energy mode, nor is simultaneous
    /// neutrino-antineutrino solution (both) possible.
    /// \details Constructors projectors and initializes Hamiltonian.
    nuSQUIDS(unsigned int numneu, NeutrinoType NT = neutrino):
    numneu(numneu),iinteraction(false),NT(NT)
    {init();}

    /// \brief Constructor from a HDF5 filepath.
    /// @param hdf5_filename Filename of the HDF5 to use for construction.
    /// @param grp HDF5 file group path where the nuSQUIDS object is saved.
    /// @param int_struct Interaction Structure to be used in object construction,
    ///        if this has already been obtained. May be left as NULL to cause
    ///        the interaction information to be read from the standard location.
    /// \details Reads the HDF5 file and constructs the associated nuSQUIDS object
    /// restoring all properties as well as the state.
    /// @see ReadStateHDF5
    nuSQUIDS(std::string hdf5_filename, std::string grp = "/",
             std::shared_ptr<InteractionStructure> int_struct = nullptr):
      int_struct(int_struct)
    {
      if(this->int_struct)
        ReadStateHDF5(hdf5_filename,grp,int_struct);
      else
        ReadStateHDF5(hdf5_filename, grp, "");
    }

    //***************************************************************
    virtual ~nuSQUIDS();

    //***************************************************************
    ///\brief Move assigns a nuSQUIDS object from an existing object
    nuSQUIDS& operator=(nuSQUIDS&&);
  
  protected:
    /************************************************************************************
     * PHYSICS FUNCTIONS - SUPER IMPORTANT
    *************************************************************************************/

    /// \brief Returns the time independent part of the Hamiltonian
    /// @param E Neutrino energy [eV]
    /// @param irho Density matrix equation index.
    /// \details \c irho is the index of the equation on which the Hamiltonian
    /// is used. When not in NeutrinoType == \c both \c irho is only 0, but if
    /// neutrinos and antineutrinos are solved simultaneously, then 0 is for
    /// neutrinos and 1 for antineutrinos.
    squids::SU_vector H0(double E, unsigned int irho) const;

    /// \brief Returns the time dependent part of the Hamiltonian at an energy node \c ie
    /// @param ie Energy node
    /// @param irho Density matrix equation index.
    /// @see H0 for convention on \c irho.
    virtual squids::SU_vector HI(unsigned int ie, unsigned int irho) const;

    /// \brief Returns noncoherent interaction term in the Hamiltonian at node \c ie
    /// @param ie Energy node
    /// @param irho Density matrix equation index.
    /// @see H0 for convention on \c irho.
    virtual squids::SU_vector GammaRho(unsigned int ie, unsigned int irho) const;

    /// \brief Returns the scalar quantity decay rate at an energy node \c ie
    /// @param ie Energy node
    /// @param iscalar scalar equation index.
    /// \details When tau regenereation is considered \c iscalar = 0 corresponds
    /// to the tau flux.
    virtual double GammaScalar(unsigned int ie, unsigned int iscalar) const;

    /// \brief Returns interaction of the density matrix at an energy node \c ie
    /// @param ie Energy node
    /// @param irho Density matrix equation index.
    /// @see H0 for convention on \c irho.
    virtual squids::SU_vector InteractionsRho(unsigned int ie, unsigned int irho) const;

    /// \brief Returns scalar interactions at an energy node \c ie
    /// @param ie Energy node
    /// @param iscalar scalar equation index.
    /// @see GammaScalar for convention on \c iscalar.
    virtual double InteractionsScalar(unsigned int ie, unsigned int iscalar) const;
  private:
    /// \brief SQuIDS signature of HI
    squids::SU_vector HI(unsigned int ie, unsigned int irho, double x) const {return HI(ie,irho);}
    /// \brief SQuIDS signature of GammaRho
    squids::SU_vector GammaRho(unsigned int ei, unsigned int irho, double x) const {return GammaRho(ei,irho);}
    /// \brief SQuIDS signature of GammaScalar
    double GammaScalar(unsigned int ei, unsigned int iscalar,double x) const {return GammaScalar(ei,iscalar);}
    /// \brief SQuIDS signature of InteractionsRho
    squids::SU_vector InteractionsRho(unsigned int ei, unsigned int irho, double x) const {return InteractionsRho(ei,irho);}
    /// \brief SQuIDS signature of InteractionsScalar
    double InteractionsScalar(unsigned int ei, unsigned int irho, double x) const {return InteractionsScalar(ei,irho);}
  public:
    /************************************************************************************
     * PUBLIC MEMBERS TO EVALUATE/SET/GET STUFF
    *************************************************************************************/

    /// \brief Sets the initial state in the single energy mode
    /// @param ini_state Initial neutrino state.
    /// @param basis Representation of the neutrino state either flavor or mass.
    /// \details \c ini_state length has to be equal to \c numneu. If the basis is 
    /// flavor then the entries are interpret as nu_e, nu_mu, nu_tau, nu_sterile_1, ..., nu_sterile_n,
    /// while if the mass basis is used then the first entries correspond to the active
    /// mass eigenstates.
    void Set_initial_state(const marray<double,1>& ini_state, Basis basis = flavor);

    /// \brief Sets the initial state in the multiple energy mode when doing either neutrino or antineutrino only.
    /// @param ini_state Initial neutrino state.
    /// @param basis Representation of the neutrino state either flavor or mass.
    /// \details \c ini_state first dimension has length equal to \c ne (number of energy nodes), while
    /// the second dimension has length equal to \c numneu (number of flavors). If the basis is
    /// flavor then the entries are interpret as nu_e, nu_mu, nu_tau, nu_sterile_1, ..., nu_sterile_n,
    /// while if the mass basis is used then the first entries correspond to the active
    /// mass eigenstates.
    void Set_initial_state(const marray<double,2>& ini_state, Basis basis = flavor);

    /// \brief Sets the initial state in the multiple energy mode when doing either neutrino or antineutrino only.
    /// @param ini_state Initial neutrino state.
    /// @param basis Representation of the neutrino state either flavor or mass.
    /// \details \c ini_state first dimension has length equal to \c ne (number of energy nodes), while
    /// the second dimension has length equal to \c numneu (number of flavors). If the basis is
    /// flavor then the entries are interpret as nu_e, nu_mu, nu_tau, nu_sterile_1, ..., nu_sterile_n,
    /// while if the mass basis is used then the first entries correspond to the active
    /// mass eigenstates.
    void Set_initial_state(const marray<double,3>& ini_state, Basis basis = flavor);

    /// \brief Sets the body where the neutrino will propagate.
    /// @param body Body object.
    void Set_Body(std::shared_ptr<Body> body);

    /// \brief Sets neutrino trajectory
    /// @param track Trajectory corresponding to the Body.
    void Set_Track(std::shared_ptr<Track> track);

    /// \brief Sets the neutrino energy in single energy mode.
    /// @param enu Neutrino energy [eV]
    /// @pre NeutrinoType must be neutrino or antineutrino, not both.
    void Set_E(double enu);

    /// \brief Evolves the system.
    /// @pre Body, track, neutrino state, and energy must has being set
    /// before calling this function.
    void EvolveState();

    /// \brief Returns the mass composition at a given node.
    /// @param flv Neutrino flavor.
    /// @param ie Energy node index.
    /// @param rho Index of the equation, see details.
    /// \details When NeutrinoType is \c both \c rho specifies wether one
    /// is considering neutrinos (0) or antineutrinos (1).

    double EvalMassAtNode(unsigned int flv,unsigned int ie,unsigned int rho = 0) const;
    /// \brief Returns the flavor composition at a given node.
    /// @param flv Neutrino flavor.
    /// @param ie Energy node index.
    /// @param rho Index of the equation, see details.
    /// \details When NeutrinoType is \c both \c rho specifies wether one
    /// is considering neutrinos (0) or antineutrinos (1).
    double EvalFlavorAtNode(unsigned int flv,unsigned int ie,unsigned int rho = 0) const;

    /// \brief Returns the mass composition at a given energy in the multiple energy mode.
    /// @param flv Neutrino flavor.
    /// @param enu Neutrino energy in natural units [eV].
    /// @param rho Index of the equation, see details.
    /// \details When NeutrinoType is \c both \c rho specifies wether one
    /// is considering neutrinos (0) or antineutrinos (1).
    double EvalMass(unsigned int flv,double enu,unsigned int rho = 0) const;

    /// \brief Returns the flavor composition at a given energy in the multiple energy mode.
    /// @param flv Neutrino flavor.
    /// @param enu Neutrino energy in natural units [eV].
    /// @param rho Index of the equation, see details.
    /// \details When NeutrinoType is \c both \c rho specifies wether one
    /// is considering neutrinos (0) or antineutrinos (1).
    double EvalFlavor(unsigned int flv,double enu,unsigned int rho = 0) const;

    /// \brief Returns the mass composition in the single energy mode.
    /// @param flv Neutrino flavor.
    double EvalMass(unsigned int flv) const;

    /// \brief Returns the flavor composition in the single energy mode.
    /// @param flv Neutrino flavor.
    double EvalFlavor(unsigned int flv) const;

    /// \brief Toggles tau regeneration on and off.
    /// \param opt If \c true tau regeneration will be considered.
    void Set_TauRegeneration(bool opt);
  
    /// \brief Toggles the effect of the Glashow resonance on and off.
    /// \param opt If \c true resonant W^- production will be considered.
    void Set_GlashowResonance(bool opt);

    /// \brief Toggles neutrino oscillations on and off.
    /// \param opt If \c true neutrino oscillations will be considered.
    void Set_IncludeOscillations(bool opt);

    /// \brief Toggles positivization of the flux.
    /// @param opt If \c true the flux will be forced to be positive every \c positivization_step.
    void Set_PositivityConstrain(bool opt);

    /// \brief Sets the positivization step.
    /// @param step Sets the positivization step.
    void Set_PositivityConstrainStep(double step);

    /// \brief Toggles the progress bar printing on and off
    /// @param opt If \c true a progress bar will be printed.
    void Set_ProgressBar(bool opt);

    /// \brief Returns the energy nodes values.
    marray<double,1> GetERange() const;

    /// \brief Returns number of energy nodes.
    size_t GetNumE() const;

    /// \brief Returns the number of neutrino flavors.
    unsigned int GetNumNeu() const;

    /// \brief Returns the number of rho equations.
    unsigned int GetNumRho() const;

    /// \brief Returns true if noncoherent interactions are considered
    bool GetUseInteractions() const {return iinteraction;}

    /// \brief Returns true if oscillations are considered
    bool GetUseOscillations() const {return ioscillations;}

    /// \brinf Initiàlized interaction structure
    void InitializeInteractions();

    /// \brief Returns the interaction structure.
    std::shared_ptr<nuSQUIDS::InteractionStructure> GetInteractionStructure() {
      if(!interactions_initialized)
        throw std::logic_error("Interactions not initialized");
      return int_struct;
    }

    /// \brief Returns the interaction structure with constness.
    std::shared_ptr<const nuSQUIDS::InteractionStructure> GetInteractionStructure() const {
      if(!interactions_initialized)
        throw std::logic_error("Interactions not initialized");
      return int_struct;
    }

    /// \brief Return the Hamiltonian at the current time
    squids::SU_vector GetHamiltonian(unsigned int ei, unsigned int rho = 0);
    /// \brief Returns the state
    const squids::SU_vector& GetState(unsigned int ei,unsigned int rho = 0) const{
      return state[ei].rho[rho];
    }
    /// \brief Returns the flavor projector
    const squids::SU_vector& GetFlavorProj(unsigned int flv,unsigned int rho = 0) const{
      return b1_proj[rho][flv];
    }
    /// \brief Returns the mass projector
    const squids::SU_vector& GetMassProj(unsigned int flv,unsigned int rho = 0) const{
      return b0_proj[flv];
    }

    /// \brief Returns the trajectory object.
    std::shared_ptr<Track> GetTrack();
    /// \brief Returns the body object.
    std::shared_ptr<Body> GetBody();

    const Track& GetTrackFast() const{ return(*track); }
  
    /// \brief Writes the object into an HDF5 file.
    /// @param hdf5_filename Filename of the HDF5 to use for construction.
    /// @param group Path to the group where the nuSQUIDS content will be saved.
    /// @param save_cross_sections If \c true the cross section tables will also be saved.
    /// @param cross_section_grp_loc Path to the group where the cross section will be saved.
    /// @param overwrite If true will overwrite previous file, else it will append.
    /// \details By default all contents are saved to the \c root of the HDF5 file
    /// if the user wants a different location it can use \c group and \c save_cross_sections to
    /// change the nuSQUIDS location and the cross section information location. Furthermore, the 
    /// \c save_cross_sections flag, which defaults to \c false decides if the cross section
    /// will be saved. Finally, it calls AddToWriteHDF5() to enable user defined parameters
    /// to be saved.
    /// @see AddToWriteHDF5
    /// @see ReadStateHDF5
    void WriteStateHDF5(std::string hdf5_filename,std::string group = "/",bool save_cross_sections = true, std::string cross_section_grp_loc = "",bool overwrite = true) const;

    /// \brief User function to write user defined properties from an HDF5 file.
    /// @param hdf5_loc_id HDF5 group id
    /// \details This function will be called by WriteStateHDF5() which will provide
    /// the hdf5_loc_id enabling the user to read user defined contents.
    /// @see WriteStateHDF5
    /// @see AddToReadHDF5
    virtual void AddToWriteHDF5(hid_t hdf5_loc_id) const;

    /// \brief Reads and constructs the object from an HDF5 file.
    /// @param hdf5_filename Filename of the HDF5 to use for construction.
    /// @param group Path to the group where the nuSQUIDS content will be saved.
    /// @param cross_section_grp_loc Path to the group where the cross section will be saved.
    /// \details By default all contents are assumed to be in the \c root of the HDF5 file
    /// if the user wants a different location it can use \c group and \c save_cross_sections to
    /// change the nuSQUIDS location and the cross section information location.
    /// Finally, it calls AddToReadHDF5() to enable user defined parameters
    /// to be read.
    /// @see AddToReadHDF5
    /// @see WriteStateHDF5
    void ReadStateHDF5(std::string hdf5_filename,std::string group = "/", std::string cross_section_grp_loc = "");


    /// \brief User function to read user defined properties from an HDF5 file.
    /// @param hdf5_loc_id HDF5 group id
    /// \details This function will be called by ReadStateHDF5() which will provide
    /// the hdf5_loc_id enabling the user to read user defined contents.
    /// @see ReadStateHDF5
    /// @see AddToWriteHDF5
    virtual void AddToReadHDF5(hid_t hdf5_loc_id);

    /// \brief Sets the mixing angle th_ij.
    /// @param i the (zero-based) index of the first state
    /// @param j the (zero-based) index of the second state
    ///              must be larger than \c i.
    /// @param angle Angle to use in radians.
    /// \details Sets the neutrino mixing angle. In our zero-based convention, e.g., the th_12 is i = 0, j = 1.,etc.
    void Set_MixingAngle(unsigned int i, unsigned int j,double angle);

    /// \brief Returns the mixing angle th_ij in radians.
    /// @param i the (zero-based) index of the first state
    /// @param j the (zero-based) index of the second state
    ///              must be larger than \c i.
    /// \details Gets the neutrino mixing angle. In our zero-based convention, e.g., the th_12 is i = 0, j = 1.,etc.
    double Get_MixingAngle(unsigned int i, unsigned int j) const;

    /// \brief Sets the CP phase for the ij-rotation.
    /// @param i the (zero-based) index of the first state
    /// @param j the (zero-based) index of the second state
    ///              must be larger than \c i.
    /// @param angle Phase to use in radians.
    /// \details Sets the CP phase for the ij-rotation. In our zero-based convention, e.g., the delta_13 = delta_CP  is i = 0, j = 2.,etc.
    void Set_CPPhase(unsigned int i, unsigned int j,double angle);

    /// \brief Returns the CP phase of the ij-rotation in radians.
    /// @param i the (zero-based) index of the first state
    /// @param j the (zero-based) index of the second state
    ///              must be larger than \c i.
    /// \details Gets the CP phase for the ij-rotation. In our zero-based convention, e.g., the delta_13 = delta_CP  is i = 0, j = 2.,etc.
    double Get_CPPhase(unsigned int i, unsigned int j) const;

    /// \brief Sets the square mass difference with respect to first mass eigenstate
    /// @param i the (zero-based) index of the second state
    ///              must be larger than \c 0.
    /// @param sq Square mass difference in eV^2.
    /// \details Sets square mass difference with respect to the first mass eigenstate. In our zero-based convention, e.g., the \f$\Delta m^2_{12}\f$ corresponds to (i = 1),etc.
    void Set_SquareMassDifference(unsigned int i, double sq);

    /// \brief Returns the square mass difference between state-i and the first mass eigenstate.
    /// @param i the (zero-based) index of the first state
    ///              must be larger than \c i.
    /// \details Returns square mass difference with respect to the first mass eigenstate. In our zero-based convention, e.g., the \f$\Delta m^2_{12}\f$ corresponds to (i = 1),etc.
    double Get_SquareMassDifference(unsigned int i) const;

    /// \brief Sets the mixing parameters to default.
    void Set_MixingParametersToDefault();

    /// \brief Sets the basis on which the evolution will be perfomed.
    /// @param basis Evolution basis can be either \c mass or \c interaction.
    /// \details By detault the solution is perfomed in the interaction basis, which
    /// is recommended for most problems. One might consider the user of the mass
    /// basis when collective effects, defined by new
    /// interactions introduce phases that cannot be cancelled as one goes to the
    /// interaction basis.
    void Set_Basis(Basis basis);

    /// \brief Sets the debug information on/off
    /// @param debug_ Boolean that toggles debuging on and off.
    void Set_Debug(bool debug_){
      debug = debug_;
    }

    /// \brief Returns lepton mixing matrix
    std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex *)> GetTransformationMatrix() const {
      return params.GetTransformationMatrix(GetNumNeu());
    }

    /// \brief Returns the neutrino interaction cross sections
    std::shared_ptr<NeutrinoCrossSections> GetNeutrinoCrossSections() {
      return(ncs);
    }
  
    /// \brief Returns the neutrino interaction cross sections
    void SetNeutrinoCrossSections(std::shared_ptr<NeutrinoCrossSections> xs) {
       ncs=xs;
       interactions_initialized=false;
    }
};

/**
 * The following class provides functionalities
 * for atmospheric neutrino experiments
 * where a collection of trayectories is explored.
 */

template<typename BaseType = nuSQUIDS, typename = typename std::enable_if<std::is_base_of<nuSQUIDS,BaseType>::value>::type >
class nuSQUIDSAtm {
  public:
    using BaseSQUIDS = BaseType;
  private:
    /// \brief Random number generator
    gsl_rng * r_gsl;
    // /// \brief Mixing matrix
    // std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex *)> U{nullptr,[](gsl_matrix_complex*m){delete m;}};
    /// \brief Internal units
    const squids::Const units;
    /// \brief Boolean that signals that a progress bar will be printed.
    bool progressbar = false;
    /// \brief Linear interpolator.
    double LinInter(double x,double xM, double xP, double yM, double yP) const {
      double f2=(x-xM)/(xP-xM);
      double f1=1-f2;
      return f1*yM + f2*yP;
    }
    /// \brief Boolean that signals that an initial state has being set.
    bool iinistate;
    /// \brief Boolean that signals the object correct initialization.
    bool inusquidsatm;
    /// \brief Contains the cos(zenith) nodes.
    marray<double,1> costh_array;
    /// \brief Contains the energy nodes [eV].
    marray<double,1> enu_array;
    /// \brief Contains the log of energy nodes (in log of eV).
    marray<double,1> log_enu_array;
    /// \brief Contains the nuSQUIDS objects for each zenith.
    std::vector<BaseSQUIDS> nusq_array;

    /// \brief Contains the Earth in atmospheric configuration.
    std::shared_ptr<EarthAtm> earth_atm;
    /// \brief Contains the trajectories for each nuSQUIDS object, i.e. zenith.
    std::vector<std::shared_ptr<EarthAtm::Track>> track_array;
    /// \brief Contains the neutrino cross section object
    std::shared_ptr<NeutrinoCrossSections> ncs;
    /// \brief Contains the interaction information structure.
    std::shared_ptr<typename BaseSQUIDS::InteractionStructure> int_struct;
  public:
    /************************************************************************************
     * CONSTRUCTORS
    *************************************************************************************/

    /// \brief Basic constructor.
    /// @param costh_array One dimensional array containing zenith angles to be calculated.
    /// @param energy_min Minimum neutrino energy value [ev].
    /// @param energy_max Maximum neutrino energy value [eV].
    /// @param energy_div Number of energy divisions.
    /// @param numneu Number of neutrino flavors.
    /// @param NT Signals the neutrino type : neutrino, antineutrion or both (simultaneous solution)
    /// @param iinteraction Sets the neutrino noncoherent neutrino interactions on.
    /// \details By defaults interactions are not considered and the neutrino energy scale is assume logarithmic.
    template<typename... ArgTypes>
    nuSQUIDSAtm(marray<double,1> costh_array, ArgTypes&&... args):
    costh_array(costh_array)
    {
      gsl_rng_env_setup();
      const gsl_rng_type * T_gsl = gsl_rng_default;
      r_gsl = gsl_rng_alloc (T_gsl);

      nusq_array = std::vector<BaseSQUIDS>(costh_array.extent(0));

      earth_atm = std::make_shared<EarthAtm>();
      for(double costh : costh_array)
        track_array.push_back(std::make_shared<EarthAtm::Track>(acos(costh)));

      unsigned int i = 0;
      for(BaseSQUIDS& nsq : nusq_array){
        nsq = BaseSQUIDS(args...);
        if(!ncs)
          ncs=nsq.GetNeutrinoCrossSections();
        nsq.Set_Body(earth_atm);
        nsq.Set_Track(track_array[i]);
        i++;
      }

      enu_array = nusq_array.front().GetERange();
      log_enu_array.resize(std::vector<size_t>{enu_array.size()});
      //std::transform(enu_array.begin(), enu_array.end(), log_enu_array.begin(),
      //               [](int enu) { return log(enu); });
      for(unsigned int ie = 0 ; ie < enu_array.size(); ie++)
        log_enu_array[ie] = log(enu_array[ie]);

      if(nusq_array.front().GetUseInteractions()){
        // setting the interaction structure also on the  main object
        nusq_array.front().InitializeInteractions();
        int_struct = nusq_array.front().GetInteractionStructure();
        for(BaseSQUIDS& nsq : nusq_array)
          nsq.int_struct=int_struct;
      }

      inusquidsatm = true;
    }

    /// \brief Constructor from a HDF5 filepath.
    /// @param hdf5_filename Filename of the HDF5 to use for construction.
    /// \details Reads the HDF5 file and construct the associated nuSQUIDSAtm object
    /// restoring all properties as well as the state.
    /// @see ReadStateHDF5
    nuSQUIDSAtm(std::string hdf5_filename) {
      gsl_rng_env_setup();
      const gsl_rng_type * T_gsl = gsl_rng_default;
      r_gsl = gsl_rng_alloc (T_gsl);
      ReadStateHDF5(hdf5_filename);
      //U = nusq_array[0].GetTransformationMatrix();
    }

    /// \brief Move constructor.
    nuSQUIDSAtm(nuSQUIDSAtm&& other):
    progressbar(other.progressbar),
    iinistate(other.iinistate),
    inusquidsatm(other.inusquidsatm),
    costh_array(std::move(other.costh_array)),
    enu_array(std::move(other.enu_array)),
    log_enu_array(std::move(other.log_enu_array)),
    nusq_array(std::move(other.nusq_array)),
    earth_atm(std::move(other.earth_atm)),
    track_array(std::move(other.track_array)),
    ncs(std::move(other.ncs)),
    int_struct(std::move(other.int_struct))
    {
      other.inusquidsatm = false;
    }

    //***************************************************************
    ///\brief Move assigns a nuSQUIDSAtm object from an existing object
    nuSQUIDSAtm& operator=(nuSQUIDSAtm&& other){
      if(&other==this)
        return(*this);

      progressbar = other.progressbar;
      iinistate = other.iinistate;
      inusquidsatm = other.inusquidsatm;
      costh_array = std::move(other.costh_array);
      enu_array = std::move(other.enu_array);
      log_enu_array = std::move(other.log_enu_array);
      nusq_array = std::move(other.nusq_array);
      earth_atm = std::move(other.earth_atm);
      track_array = std::move(other.track_array);
      ncs = std::move(other.ncs);
      int_struct = std::move(other.int_struct);

      // initial nusquids object render useless
      other.inusquidsatm = false;

      return(*this);
    }
    /************************************************************************************
     * PUBLIC MEMBERS TO EVALUATE/SET/GET STUFF
    *************************************************************************************/

    /// \brief Sets the initial state in the multiple energy mode
    /// when only considering neutrinos or antineutrinos
    /// @param ini_state Initial neutrino state.
    /// @param basis Representation of the neutrino state either flavor or mass.
    /// \details \c ini_state first dimension lenght has to be equal to the number of
    /// zenith bins, the second dimension is the number of energy bins, and the third
    /// dimension corresponds to the number of neutrino flavors. If the basis is
    /// flavor then the entries are interpret as nu_e, nu_mu, nu_tau, nu_sterile_1, ..., nu_sterile_n,
    /// while if the mass basis is used then the first entries correspond to the active
    /// mass eigenstates.
    void Set_initial_state(marray<double,3> ini_flux, Basis basis){
      if(ini_flux.extent(0) != costh_array.extent(0))
        throw std::runtime_error("nuSQUIDSAtm::Error::First dimension of input array is incorrect.");
      if(ini_flux.extent(1) != enu_array.extent(0))
        throw std::runtime_error("nuSQUIDSAtm::Error::Second dimension of input array is incorrect.");
      unsigned int i = 0;
      for(BaseSQUIDS& nsq : nusq_array){
        marray<double,2> slice{ini_flux.extent(1),ini_flux.extent(2)};
        for(size_t j=0; j<ini_flux.extent(1); j++){
          for(size_t k=0; k<ini_flux.extent(2); k++)
            slice[j][k]=ini_flux[i][j][k];
        }
        nsq.Set_initial_state(slice,basis);
        //nsq.Set_initial_state(ini_flux[i],basis);
        i++;
      }
      iinistate = true;
    }

    /// \brief Sets the initial state in the multiple energy mode when
    /// considering neutrinos and antineutrinos simultaneously
    /// @param ini_state Initial neutrino state.
    /// @param basis Representation of the neutrino state either flavor or mass.
    /// \details \c ini_state first dimension lenght has to be equal to the number of
    /// zenith bins, the second dimension is the number of energy bins, the third
    /// dimension must have length two corresponding to neutrinos and antineutrinos 
    /// respectifully, finally the fourth
    /// dimension corresponds to the number of neutrino flavors. If the basis is
    /// flavor then the entries are interpret as nu_e, nu_mu, nu_tau, nu_sterile_1, ..., nu_sterile_n,
    /// while if the mass basis is used then the first entries correspond to the active
    /// mass eigenstates.
    void Set_initial_state(marray<double,4> ini_flux, Basis basis){
      if(ini_flux.extent(0) != costh_array.extent(0))
        throw std::runtime_error(
            "nuSQUIDSAtm::Error::First dimension of input array is incorrect.");
      unsigned int i = 0;
      for(BaseSQUIDS& nsq : nusq_array){
        marray<double,3> slice{ini_flux.extent(1),ini_flux.extent(2),ini_flux.extent(3)};
        for(size_t j=0; j<ini_flux.extent(1); j++){
          for(size_t k=0; k<ini_flux.extent(2); k++){
            for(size_t m=0; m<ini_flux.extent(3); m++)
            slice[j][k][m]=ini_flux[i][j][k][m];
          }
        }

        nsq.Set_initial_state(slice,basis);
        //nsq.Set_initial_state(ini_flux[i],basis);
        i++;
      }
      iinistate = true;
    }

    /// \brief Evolves the system.
    void EvolveState(){
      if(not iinistate)
        throw std::runtime_error("nuSQUIDSAtm::Error::State not initialized.");
      if(not inusquidsatm)
        throw std::runtime_error("nuSQUIDSAtm::Error::nuSQUIDSAtm not initialized.");
      unsigned int i = 0;
      for(BaseSQUIDS& nsq : nusq_array){
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

    /// \brief Returns the flavor composition at a given energy and zenith.
    /// @param flv Neutrino flavor.
    /// @param costh Cosine of the zenith.
    /// @param enu Neutrino energy [eV].
    /// @param rho Index of the equation, see details.
    /// \details When NeutrinoType is \c both \c rho specifies wether one
    /// is considering neutrinos (0) or antineutrinos (1). Bilinear interpolation
    /// is done in the logarithm of the energy and cos(zenith).
    double EvalFlavor(unsigned int flv,double costh,double enu,unsigned int rho = 0, bool randomize_production_height = false) const {
      // here the energy enters in eV
      if(not iinistate)
        throw std::runtime_error("nuSQUIDSAtm::Error::State not initialized.");
      if(not inusquidsatm)
        throw std::runtime_error("nuSQUIDSAtm::Error::nuSQUIDSAtm not initialized.");
      
      if( costh < *costh_array.begin() or costh > *costh_array.rbegin())
        throw std::runtime_error("nuSQUIDSAtm::Error::EvalFlavor::cos(th) out of bounds.");
      if( enu < *enu_array.begin() or enu > *enu_array.rbegin() )
        throw std::runtime_error("nuSQUIDSAtm::Error::EvalFlavor::neutrino energy out of bounds.(Emin = " +
                                 std::to_string(*enu_array.begin()) +
                                 ",Emax = " +
                                 std::to_string(*enu_array.rbegin()) +
                                 ", Enu = " + std::to_string(enu) + ")");
      
      auto cthit=std::lower_bound(costh_array.begin(),costh_array.end(),costh);
      if(cthit==costh_array.end())
        throw std::runtime_error("SQUIDS::GetExpectationValueD : x value not in the array.");
      if(cthit!=costh_array.begin())
        cthit--;
      size_t cth_M=std::distance(costh_array.begin(),cthit);
      
      /*double logE = log(enu);
      auto logeit=std::lower_bound(log_enu_array.begin(),log_enu_array.end(),logE);
      if(logeit==log_enu_array.end())
        throw std::runtime_error("SQUIDS::GetExpectationValueD : x value not in the array.");
      if(logeit!=log_enu_array.begin())
        logeit--;
      size_t loge_M=std::distance(log_enu_array.begin(),logeit);*/
      
      auto eit=std::lower_bound(enu_array.begin(),enu_array.end(),enu);
      if(eit==enu_array.end())
        throw std::runtime_error("SQUIDS::GetExpectationValueD : x value not in the array.");
      if(eit!=enu_array.begin())
        eit--;
      size_t loge_M=std::distance(enu_array.begin(),eit);
      
      //EarthAtm::Track track(acos(costh));
      EarthAtm::Track track=EarthAtm::Track::makeWithCosine(costh);
      double delta_t_final = track.GetFinalX()-track.GetInitialX();
      if (randomize_production_height){
        double production_height = gsl_ran_flat(r_gsl,-15*units.km,15*units.km);
        delta_t_final += production_height;
      }
      
      // assuming offsets are zero
      double delta_t_1 = nusq_array[cth_M].Get_t() - nusq_array[cth_M].Get_t_initial();
      double delta_t_2 = nusq_array[cth_M+1].Get_t() - nusq_array[cth_M+1].Get_t_initial();
      double delta_t_final_1 = nusq_array[cth_M].GetTrackFast().GetFinalX() - nusq_array[cth_M].GetTrackFast().GetInitialX();
      double delta_t_final_2 = nusq_array[cth_M+1].GetTrackFast().GetFinalX() - nusq_array[cth_M+1].GetTrackFast().GetInitialX();
      double t_inter = 0.5*(delta_t_final*delta_t_1/delta_t_final_1 + delta_t_final*delta_t_2/delta_t_final_2);
      
      squids::SU_vector H0_at_enu = nusq_array[0].H0(enu,rho);
      
      struct storage_type{
        squids::SU_vector evol_proj, temp1, temp2;
        storage_type(unsigned int dim):evol_proj(dim),temp1(dim),temp2(dim){}
      };
  #ifdef SQUIDS_THREAD_LOCAL
      static SQUIDS_THREAD_LOCAL
  #endif
      storage_type storage(H0_at_enu.Dim());
      
      storage.evol_proj = nusq_array[0].GetFlavorProj(flv,rho).Evolve(H0_at_enu,t_inter);
      
      //coefficients for energy interpolation
      double f2=(enu-enu_array[loge_M])/(enu_array[loge_M+1]-enu_array[loge_M]);
      double f1=1-f2;
      
      storage.temp1 =f1*nusq_array[cth_M  ].GetState(loge_M  ,rho);
      storage.temp1+=f2*nusq_array[cth_M  ].GetState(loge_M+1,rho);
      storage.temp2=storage.temp1.Evolve(H0_at_enu,t_inter - nusq_array[cth_M  ].Get_t());
      double phiM=squids::SUTrace<squids::detail::AlignedStorage>(storage.temp2,storage.evol_proj);
      
      storage.temp1 =f1*nusq_array[cth_M+1].GetState(loge_M  ,rho);
      storage.temp1+=f2*nusq_array[cth_M+1].GetState(loge_M+1,rho);
      storage.temp2=storage.temp1.Evolve(H0_at_enu,t_inter - nusq_array[cth_M+1].Get_t());
      double phiP=squids::SUTrace<squids::detail::AlignedStorage>(storage.temp2,storage.evol_proj);
      
      //perform angular interpolation
      return(LinInter(costh,costh_array[cth_M],costh_array[cth_M+1],phiM,phiP));
    }

    /// \brief Writes the object into an HDF5 file.
    /// @param hdf5_filename Filename of the HDF5 into which save the object.
    /// \details All contents are saved to the \c root of the HDF5 file.
    /// @see ReadStateHDF5
    void WriteStateHDF5(std::string filename, bool overwrite = true) const{
      if(not iinistate)
        throw std::runtime_error("nuSQUIDSAtm::Error::State not initialized.");
      if(not inusquidsatm)
        throw std::runtime_error("nuSQUIDSAtm::Error::nuSQUIDSAtm not initialized.");

      hid_t file_id,root_id;
      hid_t dset_id;
      // create HDF5 file
      if(overwrite)
        file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      else
        file_id = H5Fcreate (filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT, H5P_DEFAULT);
      if (file_id < 0)
          throw std::runtime_error("nuSQUIDS::Error::Cannot create file at " + filename + ".");
      root_id = H5Gopen(file_id, "/",H5P_DEFAULT);

      // write the zenith range
      hsize_t costhdims[1]={costh_array.extent(0)};
      dset_id = H5LTmake_dataset(root_id,"zenith_angles",1,costhdims,H5T_NATIVE_DOUBLE,costh_array.get_data());
      hsize_t energydims[1]={enu_array.extent(0)};
      dset_id = H5LTmake_dataset(root_id,"energy_range",1,energydims,H5T_NATIVE_DOUBLE,enu_array.get_data());

      H5Gclose (root_id);
      H5Fclose (file_id);

      unsigned int i = 0;
      for(const BaseSQUIDS& nsq : nusq_array){
        // use only the first one to write the cross sections on /crosssections
        // the last false is because the HDF5 file has already been created, should not overwrite 
        nsq.WriteStateHDF5(filename,"/costh_"+std::to_string(costh_array[i]),i==0,"crosssections",false);
        i++;
      }
    }

    /// \brief Reads the object from an HDF5 file.
    /// @param hdf5_filename Filename of the HDF5 to use for construction.
    /// \details All contents are assumed to be saved to the \c root of the HDF5 file.
    /// @see WriteStateHDF5
    void ReadStateHDF5(std::string hdf5_filename){
      hid_t file_id,group_id,root_id;
      // create HDF5 file
      file_id = H5Fopen(hdf5_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if (file_id < 0)
        throw std::runtime_error("nuSQUIDSAtm::ReadStateHDF5: Unable to open file: " + hdf5_filename + ". No such file or directory");
      root_id = H5Gopen(file_id, "/",H5P_DEFAULT);
      group_id = root_id;

      // read the zenith range dimension
      hsize_t costhdims[1];
      H5LTget_dataset_info(group_id, "zenith_angles", costhdims, NULL, NULL);

      double data[costhdims[0]];
      H5LTread_dataset_double(group_id, "zenith_angles", data);
      costh_array.resize(std::vector<size_t> {costhdims[0]});
      for (unsigned int i = 0; i < costhdims[0]; i ++)
        costh_array[i] = data[i];

      hsize_t energydims[1];
      H5LTget_dataset_info(group_id, "energy_range", energydims, NULL, NULL);

      double enu_data[energydims[0]];
      H5LTread_dataset_double(group_id, "energy_range", enu_data);
      enu_array.resize(std::vector<size_t>{energydims[0]});log_enu_array.resize(std::vector<size_t>{energydims[0]});
      for (unsigned int i = 0; i < energydims[0]; i ++){
        enu_array[i] = enu_data[i];
        log_enu_array[i] = log(enu_data[i]);
      }

      H5Gclose(root_id);
      H5Fclose(file_id);

      // resize apropiately the nuSQUIDSAtm container vector
      nusq_array.clear();
      nusq_array = std::vector<BaseSQUIDS>(costhdims[0]);

      unsigned int i = 0;
      for(BaseSQUIDS& nsq : nusq_array){
        if(i==0){
          // read the cross sections stored in /crosssections
          nsq.ReadStateHDF5(hdf5_filename,"/costh_"+std::to_string(costh_array[i]),"crosssections");
          if(nsq.iinteraction)
            int_struct = nsq.GetInteractionStructure();
        } else {
          // read the cross sections stored in /crosssections
          nsq.ReadStateHDF5(hdf5_filename,"/costh_"+std::to_string(costh_array[i]),int_struct);
        }
        i++;
      }

      iinistate = true;
      inusquidsatm = true;
    }

    /// \brief Sets the mixing parameters to default.
    void Set_MixingParametersToDefault(){
      for(BaseSQUIDS& nsq : nusq_array){
        nsq.Set_MixingParametersToDefault();
      }
    }

    /// \brief Sets the mixing angle th_ij.
    /// @param i the (zero-based) index of the first state
    /// @param j the (zero-based) index of the second state
    ///              must be larger than \c i.
    /// @param angle Angle to use in radians.
    /// \details Sets the neutrino mixing angle. In our zero-based convention, e.g., the th_12 is i = 0, j = 1.,etc.
    void Set_MixingAngle(unsigned int i, unsigned int j,double angle){
      for(BaseSQUIDS& nsq : nusq_array){
        nsq.Set_MixingAngle(i,j,angle);
      }
    }

    /// \brief Returns the mixing angle th_ij in radians.
    /// @param i the (zero-based) index of the first state
    /// @param j the (zero-based) index of the second state
    ///              must be larger than \c i.
    /// \details Gets the neutrino mixing angle. In our zero-based convention, e.g., the th_12 is i = 0, j = 1.,etc.
    double Get_MixingAngle(unsigned int i, unsigned int j) const{
      return nusq_array[0].Get_MixingAngle(i,j);
    }

    /// \brief Sets the CP phase for the ij-rotation.
    /// @param i the (zero-based) index of the first state
    /// @param j the (zero-based) index of the second state
    ///              must be larger than \c i.
    /// @param angle Phase to use in radians.
    /// \details Sets the CP phase for the ij-rotation. In our zero-based convention, e.g., the delta_13 = delta_CP  is i = 0, j = 2.,etc.
    void Set_CPPhase(unsigned int i, unsigned int j,double angle){
      for(BaseSQUIDS& nsq : nusq_array){
        nsq.Set_CPPhase(i,j,angle);
      }
    }

    /// \brief Returns the CP phase of the ij-rotation in radians.
    /// @param i the (zero-based) index of the first state
    /// @param j the (zero-based) index of the second state
    ///              must be larger than \c i.
    /// \details Gets the CP phase for the ij-rotation. In our zero-based convention, e.g., the delta_13 = delta_CP  is i = 0, j = 2.,etc.
    double Get_CPPhase(unsigned int i, unsigned int j) const{
      return nusq_array[0].Get_CPPhase(i,j);
    }

    /// \brief Sets the square mass difference with respect to first mass eigenstate
    /// @param i the (zero-based) index of the second state
    ///              must be larger than \c 0.
    /// @param sq Square mass difference in eV^2.
    /// \details Sets square mass difference with respect to the first mass eigenstate. In our zero-based convention, e.g., the \f$\Delta m^2_{12}\f$ corresponds to (i = 1),etc.
    void Set_SquareMassDifference(unsigned int i,double sq){
      for(BaseSQUIDS& nsq : nusq_array){
        nsq.Set_SquareMassDifference(i,sq);
      }
    }

    /// \brief Returns the square mass difference between state-i and the first mass eigenstate.
    /// @param i the (zero-based) index of the first state
    ///              must be larger than \c i.
    /// \details Returns square mass difference with respect to the first mass eigenstate. In our zero-based convention, e.g., the \f$\Delta m^2_{12}\f$ corresponds to (i = 1),etc.
    double Get_SquareMassDifference(unsigned int i) const{
      return nusq_array[0].Get_SquareMassDifference(i);
    }

    /// \brief Set the initial integration step size
    /// \param h initial step size
    void Set_h(double h){
      for(BaseSQUIDS& nsq : nusq_array)
        nsq.Set_h(h);
    }
    
    /// \brief Set the initial integration step size for a single energy node
    /// \param h initial step size
    /// \param idx energy node index
    void Set_h(double h, unsigned int idx){
      nusq_array[idx].Set_h(h);
    }
  
    /// \brief Set the maximum integration step size
    /// \param h maximum step size
    void Set_h_max(double h){
      for(BaseSQUIDS& nsq : nusq_array)
        nsq.Set_h_max(h);
    }
  
    /// \brief Set the maximum integration step size for a single energy node
    /// \param h maximum step size
    /// \param idx energy node index
    void Set_h_max(double h, unsigned int idx){
      nusq_array[idx].Set_h_max(h);
    }
  
    /// \brief Set the minimum integration step size
    /// \param h minimum step size
    void Set_h_min(double h){
      for(BaseSQUIDS& nsq : nusq_array)
        nsq.Set_h_min(h);
    }
    
    /// \brief Set the minimum integration step size for a single energy node
    /// \param h minimum step size
    /// \param idx energy node index
    void Set_h_min(double h, unsigned int idx){
      nusq_array[idx].Set_h_min(h);
    }

    /// \brief Set the allowed absolute numerical error
    /// \param eps maximum allowed error
    void Set_abs_error(double eps){
      for(BaseSQUIDS& nsq : nusq_array)
        nsq.Set_abs_error(eps);
    }
  
    /// \brief Set the allowed absolute numerical error for a single energy node
    /// \param eps maximum allowed error
    /// \param idx energy node index
    void Set_abs_error(double eps, unsigned int idx){
      nusq_array[idx].Set_abs_error(eps);
    }

    /// \brief Set the allowed relative numerical error
    /// \param eps maximum allowed error
    void Set_rel_error(double eps){
      for(BaseSQUIDS& nsq : nusq_array)
        nsq.Set_rel_error(eps);
    }
  
    /// \brief Set the allowed relative numerical error for a single energy node
    /// \param eps maximum allowed error
    /// \param idx energy node index
    void Set_rel_error(double eps, unsigned int idx){
      nusq_array[idx].Set_rel_error(eps);
    }

    /// \brief Sets the GSL solver
    /// @param opt GSL stepper function.
    void Set_GSL_step(gsl_odeiv2_step_type const * opt){
      for(BaseSQUIDS& nsq : nusq_array){
        nsq.Set_GSL_step(opt);
      }
    }

    /// \brief Toggles the progress bar printing on and off
    /// @param opt If \c true a progress bar will be printed.
    void Set_ProgressBar(bool opt){
        progressbar = opt;
        for(BaseSQUIDS& nsq : nusq_array){
          nsq.Set_ProgressBar(opt);
        }
    }

    /// \brief Returns the interaction structure.
    std::shared_ptr<const nuSQUIDS::InteractionStructure> GetInteractionStructure() const {
      return int_struct;
    }

    /// \brief Returns number of energy nodes.
    size_t GetNumE() const{
      return enu_array.extent(0);
    }

    /// \brief Returns number of zenith nodes.
    size_t GetNumCos() const{
      return costh_array.extent(0);
    }

    /// \brief Returns the number of neutrino flavors.
    unsigned int GetNumNeu() const{
      return nusq_array[0].GetNumNeu();
    }

    /// \brief Returns the number of rho equations.
    unsigned int GetNumRho() const{
      return nusq_array[0].GetNumRho();
    }

    /// \brief Returns the energy nodes values.
    marray<double,1> GetERange() const{
      return enu_array;
    }

    /// \brief Returns the cos(zenith) nodes values.
    marray<double,1> GetCosthRange() const{
      return costh_array;
    }

    /// \brief Returns the nuSQUIDSBase object for ci-th zenith.
    BaseSQUIDS& GetnuSQuIDS(unsigned int ci) {
      if(ci >= nusq_array.size())
        throw std::runtime_error("nuSQUIDSAtm::GetnuSQuIDS: input index out of range.");
      return nusq_array[ci];
    }

    /// \brief Returns the vector of nuSQUIDSBase objects.
    std::vector<BaseSQUIDS>& GetnuSQuIDS() {
      return nusq_array;
    }

    /// \brief Toggles tau regeneration on and off.
    /// @param opt If \c true tau regeneration will be considered.
    void Set_TauRegeneration(bool opt){
      for(BaseSQUIDS& nsq : nusq_array){
        nsq.Set_TauRegeneration(opt);
      }
    }

    /// \brief Toggles neutrino oscillations on and off.
    /// @param opt If \c true neutrino oscillations will be considered.
    void Set_IncludeOscillations(bool opt){
      for(BaseSQUIDS& nsq : nusq_array){
        nsq.Set_IncludeOscillations(opt);
      }
    }

    /// \brief Toggles positivization of the flux.
    /// @param opt If \c true the flux will be forced to be positive every \c positivization_step.
    void Set_PositivityConstrain(bool opt){
      for(BaseSQUIDS& nsq : nusq_array){
        nsq.Set_PositivityConstrain(opt);
      }
    }
    /// \brief Stes the step upon which the positivity correction would be apply.
    /// @param step The step upon which the positivization will take place.
    void Set_PositivityConstrainStep(double step){
      for(BaseSQUIDS& nsq : nusq_array){
        nsq.Set_PositivityConstrainStep(step);
      }
    }
};


} // close namespace
#endif
