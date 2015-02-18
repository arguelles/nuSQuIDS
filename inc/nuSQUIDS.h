#ifndef __nuSQUID_H
#define __nuSQUID_H

#include "body.h"
#include "xsections.h"
#include "taudecay.h"
#include "marray.h"

#include <algorithm>
#include <SQuIDS/SQUIDS.h>
#include <memory>
#include <map>
#include <stdexcept>

#include "H5Epublic.h"
#include "H5Tpublic.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Gpublic.h"
#include "H5Fpublic.h"

#define ManualTauReinjection
#define FixCrossSections

//#define UpdateInteractions_DEBUG

namespace nusquids{

enum MixingParameter {
  TH12 = 0,
  TH13 = 1, TH23 = 2,
  TH14 = 3, TH24 = 4, TH34 = 5,
  TH15 = 6, TH25 = 7, TH35 = 8, TH45 = 9,
  TH16 = 10, TH26 = 11, TH36 = 12, TH46 = 13, TH56 = 14,
  DELTA1 = 15, DELTA2 = 16, DELTA3 = 17,
  DM21SQ = 18, DM31SQ = 19, DM41SQ = 20, DM51SQ = 21, DM61SQ = 22
                    };

enum BASIS { mass, interaction };

static std::map<int,std::string> param_label_map {
  {0 , "th12"},
  {1 , "th13"}, {2 , "th23"},
  {3 , "th14"}, {4 , "th24"}, {5 , "th34"},
  {6 , "th15"}, {7 , "th25"}, {8 , "th35"}, {9 , "th45"},
  {10 , "th16"}, {11 , "th26"}, {12 , "th36"}, {13 , "th46"}, {14 , "th56"},
  {15 , "delta1"} , {16, "delta2"} , {17, "delta3"},
  {18, "dm21sq"} , {19 , "dm31sq"}, {20,"dm41sq"} ,{ 21 , "dm51sq"} , {22, "dm61sq"}
                                     };

static std::map<int,std::vector<int>> param_label_index {
  {0 , {0,1}},
  {1 , {0,2}}, {2 , {1,2}},
  {3 , {0,3}}, {4 , {1,3}}, {5 , {2,3}},
  {6 , {0,4}}, {7 , {1,4}}, {8 , {2,4}}, {9 , {3,4}},
  {10 , {0,5}},{11 , {1,5}}, {12 , {2,5}}, {13 , {3,5}}, {14 , {4,5}},
  {15 , {0,2}} , {16, {0,3}} , {17, {0,4}},
  {18, {1,0}} , {19 , {2,0}}, {20,{3,0}} ,{ 21 , {4,0}} , {22, {5,0}}
                                                    };
enum NeutrinoType {
  neutrino=0b01,
  antineutrino=0b10,
  both=0b11
};

///\brief nu-SQuIDS main class
class nuSQUIDS: public SQUIDS {
  // nuSQUIDSAtm is a friend class so it can use H0 to evolve operators
  // and thus evaluate expectation values.
  friend class nuSQUIDSAtm;
  protected:
    /// \brief Sets the basis in which the problem will be solved.
    ///
    /// If interaction basis is used the projectors will be evolved at
    /// every time step. On the other hand, if mass basis is used no evolution
    /// is performed.
    BASIS basis = interaction;
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
    NeutrinoCrossSections ncs;
    /// \brief Neutrino charge current differential cross section with respect to
    /// the outgoing lepton energy.
    ///
    /// The array four dimensional is constructed when InitializeInteractionVectors() is called and
    /// its initialized when InitializeInteractions() is called. The first dimension
    /// is number of neutrino types (neutrino/antineutrino/both), the second the neutrino flavor,
    /// and the last two the initial and final energy node respectively.
    marray<double,4> dNdE_CC;
    /// \brief Neutrino neutral current differential cross section with respect to
    /// the outgoing lepton energy.
    ///
    /// The array four dimensional is constructed when InitializeInteractionVectors() is called and
    /// its initialized when InitializeInteractions() is called. The first dimension
    /// is number of neutrino types (neutrino/antineutrino/both), the second the neutrino flavor,
    /// and the last two the initial and final energy node respectively.
    marray<double,4> dNdE_NC;
    /// \brief Array that contains the inverse of the neutrino neutral current mean free path.
    /// \details The array contents are in natural units (i.e. eV) and is update when
    /// UpdateInteractions() is called. The first dimension corresponds to the neutrino type,
    /// the second to the flavor, and the last one to the energy.
    marray<double,3> invlen_NC;
    /// \brief Array that contains the inverse of the neutrino charge current mean free path.
    /// \details The array contents are in natural units (i.e. eV) and is update when
    /// UpdateInteractions() is called. The first dimension corresponds to the neutrino type,
    /// the second to the flavor, and the last one to the energy.
    marray<double,3> invlen_CC;
    /// \brief Array that contains the inverse of the neutrino total mean free path.
    /// \details The array contents are in natural units (i.e. eV) and is update when
    /// UpdateInteractions() is called. Numerically it is just nuSQUIDS::invlen_NC and nuSQUIDS::invlen_CC
    /// added together.
    marray<double,3> invlen_INT;
    /// \brief Array that contains the neutrino charge current cross section.
    /// \details The first dimension corresponds to the neutrino type, the second to the flavor, and
    /// the final one to the energy node. Its contents are in natural units, i.e. eV^-2. It is
    /// initialized by InitializeInteractions() .
    marray<double,3> sigma_CC;
    /// \brief Array that contains the neutrino neutral current cross section.
    /// \details The first dimension corresponds to the neutrino type, the second to the flavor, and
    /// the final one to the energy node. Its contents are in natural units, i.e. eV^-2. It is
    /// initialized by InitializeInteractions() .
    marray<double,3> sigma_NC;

    /// \brief Interface that calculate and interpolates tau decay spectral functions.
    TauDecaySpectra tdc;
    marray<double,1> invlen_tau;
    marray<double,2> dNdE_tau_all;
    marray<double,2> dNdE_tau_lep;
    /// \brief Tau branching ratio to leptons.
    double taubr_lep;
    /// \brief Tau lifetime in natural units.
    double tau_lifetime;
    /// \brief Tau mass in natural units.
    double tau_mass;
    /// \brief Length upon which charge tau lepton conversion to neutrino happens.
    /// \defailts By default set to 100 km (in natural units).
    /// @see ConvertTauIntoNuTau()
    double tau_reg_scale;

    /// \brief Body where the neutrino propagation takes place.
    std::shared_ptr<Body> body;
    /// \brief Trayectory within the body.
    /// \details Stores the position within the body and its updated every evolution
    /// step.
    std::shared_ptr<Track> track;

    /// \brief SU_vector that represents the neutrino square mass difference matrix in the mass basis.
    ///  It is used to construct nuSQUIDS#H0_array and H0()
    SU_vector DM2;
    /// \brief Stores the time independent hamiltonian corresponding to each energy node.
    marray<SU_vector,1> H0_array;

    /// \brief Mass basis projectors.
    /// \details The i-entry corresponds to the projector in the ith mass eigenstate.
    marray<SU_vector,1> b0_proj;
    /// \brief Flavor basis projectors.
    /// \details The first dismension corresponds to the neutrino type. When NeutrinoType = both, then
    /// the first dimension equal 0 means neutrinos and 1 antineutrinos. The second dimension corresponds
    /// to the flavor eigenstates where 0 corresponds to nu_e, 1 to nu_mu, 2 to nu_tau, and the others
    /// to sterile neutrinos.
    marray<SU_vector,2> b1_proj;
    /// \brief Mass basis projectors in the interaction picture.
    /// The index meaning are the same as nuSQUIDS#b1_proj but for mass eigenstates,
    /// with an added third dimension that corresponds to the energy node index.
    marray<SU_vector,3> evol_b0_proj;
    /// \brief Flavor basis projectors in the interaction picture.
    /// The index meaning are the same as nuSQUIDS#b1_proj ,
    /// with an added third dimension that corresponds to the energy node index.
    marray<SU_vector,3> evol_b1_proj;

    /// \brief Evolves the flavor projectors in the interaction basis to a time t.
    /// \details It uses H0() to evolve SQUIDS#b0_proj and SQUIDS#b1_proj into 
    /// SQUIDS#evol_b0_proj and SQUIDS#evol_b1_proj.
    /// \warning Since the RHS of the differential equation only involves flavor projectors
    /// we do not current evolve mass projectors.
    void EvolveProjectors(double t);

    /// \brief When called converts the flux of tau charged leptons to neutrinos. It
    /// is only called when tauregeneration = True and NeutrinoType = both.
    ///
    /// At a given time we the system has a flux of charged tau leptons \phi_\tau(E) and
    /// \phi_\taubar(E) where E is given in the system nodes. When the tau/taubar decays
    /// it will produce a taunu/taunubar following a distribution dN/dE_all, where all
    /// means that both hadronic and leptonic decay channels are considered.
    /// Furthermore, the taus can decay leptonically with a distribution dN/dE_lep
    /// and branching ration br_lep.
    ///
    /// Considering this we calculate the resulting neutrinos fluxes and add them
    /// to the current neutrino state flux. In order to do this it is useful to define
    /// the following functions:
    ///
    /// \f$ F_{\tau|\bar{\tau}}(E) = \int_E^\infty d\bar{E} \frac{dN}{d\bar{E}_{all}} \phi_{\tau|\bar{\tau}}\f$
    ///
    /// then
    ///
    /// \f{eqnarray*}{
    /// \phi_{\nu_\tau}(E) &=& F_{\tau}\, \\
    /// \phi_{\bar{\nu}_\tau}(E) &=& F_{\bar{\tau}}\,\\
    /// \phi_{\nu_e} &=& F_{\bar{\tau}}*br_{lep}\,\\
    /// \phi_{\bar{\nu}_e} &=& F_{\tau}*br_{lep},
    /// \f}
    ///
    /// and similarly for \f$\nu_\mu\f$. After this transformation is done SetScalarsToZero()
    /// is called to remove the charged lepton flux.
    void ConvertTauIntoNuTau();

    // bool requirements
  private:
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
    bool iinteraction = false;
    /// \brief When multienergy mode is used, it signals that the neutrino energies is in logarithmic scale.
    bool elogscale = true;
    /// \brief Boolean that signals that tau regeneration is being used.
    bool tauregeneration = false;
    /// \brief Boolean that signals that a progress bar will be printed.
    bool progressbar = false;
    /// \brief Integer to keep track of the progress bar evolution.
    int progressbar_count = 0;
    /// \brief Number of steps upon which the progress bar will be updated.
    int progressbar_loop = 100;
    /// \brief NT keeps track if the problem consists of neutrinos, antineutrinos, or both.
    NeutrinoType NT = both;
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
    /// \defails It constructs nuSQUIDS#DM2 and nuSQUIDS#H0_array
    void iniH0();
    /// \brief Changes the CP sign of the complex matrices stored in params.
    void AntineutrinoCPFix(unsigned int irho);
    /// \brief Serializes the initialization of the Body and Track objects.
    /// @see ReadStateHDF5
    void SetBodyTrack(int,int,double*,int,double*);

    /// \brief General initilizer for the multi energy mode
    /// @param Emin Minimum neutrino energy [GeV].
    /// @param Emax Maximum neutirno energy [GeV].
    /// @param Esize Number of energy nodes.
    /// @param initialize_intereractions Togles interaction arrays initialization.
    /// @param xini The initial position of the system.
    /// \defails Constructs the energy node arrays from that energy range; the variable nuSQUIDS#elogscale
    /// sets if the scale will be linear or logarithmic. If \c initialize_intereractions 
    /// is set to \c true then InitializeInteractionVectors() and InitializeInteractions()
    /// are called.
    void init(double Emin,double Emax,unsigned int Esize,bool initialize_intereractions = true, double xini = 0.0);
    /// \brief Initilizer for the single energy mode
    /// @param xini The initial position of the system. By default is set to 0.
    void init(double xini = 0.0);
    /// \brief Initilizes auxiliary cross section arrays.
    /// \details The arrays are initialize, but not filled with contented. To fill the arrays
    /// call InitializeInteractions().
    void InitializeInteractionVectors();
    /// \brief Fills in auxiliary cross section arrays.
    /// \details It uses nuSQUIDS#ncs and nuSQUIDS#tdc to fill in the values of
    /// nuSQUIDS#dNdE_CC , nuSQUIDS#dNdE_NC , nuSQUIDS#sigma_CC , nuSQUIDS#sigma_NC ,
    /// nuSQUIDS#invlen_CC , nuSQUIDS#invlen_NC , nuSQUIDS#invlen_INT ,
    /// nuSQUIDS#dNdE_tau_all , nuSQUIDS#dNdE_tau_lep , nuSQUIDS#invlen_tau
    /// @see InitializeInteractionVectors
    void InitializeInteractions();
  private:

    /// \brief Sets all scalar arrays to zero.
    void SetScalarsToZero();
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
    virtual void AddToPreDerive(double x){};
  public:
    /// \brief Incorporated const object useful to evaluate units.
    const Const units;
    /************************************************************************************
     * CONSTRUCTORS
    *************************************************************************************/

    /// \brief Default void constructor.
    nuSQUIDS(){};
    /// \brief Multiple energy mode constructor.
    /// @param Emin Minimum neutrino energy [GeV].
    /// @param Emax Maximum neutirno energy [GeV].
    /// @param Esize Number of energy nodes.
    /// @param numneu Number of neutrino flavors.
    /// @param NT NeutrinoType: neutrino,antineutrino, or both (simultaneous solution).
    /// @param elogscale Sets the energy scale to be logarithmic
    /// @param iinteraction Sets the neutrino noncoherent neutrino interactions on.
    /// \details By default the energy scale is logarithmic and interactions are turn off.
    /// \warning When interactions are present interpolation is performed to precalculate the neutrino
    /// cross section which make take considertable time depending on the energy grid.
    /// @see init
    nuSQUIDS(double Emin,double Emax,unsigned int Esize,unsigned int numneu,NeutrinoType NT = both,
       bool elogscale = true,bool iinteraction = false):
    iinteraction(iinteraction),elogscale(elogscale),numneu(numneu),NT(NT)
    {init(Emin,Emax,Esize);};

    /// \brief Single energy mode constructor.
    /// @param numneu Number of neutrino flavors.
    /// @param NT NeutrinoType: neutrino or antineutrino.
    /// \warning Interactions are not possible in the single energy mode, nor is simultaneous
    /// neutrino-antineutrino solution (both) possible.
    /// \details Constructors projectors and initializes Hamiltonian.
    nuSQUIDS(int numneu, NeutrinoType NT = neutrino):
    iinteraction(false),elogscale(false),numneu(numneu),NT(NT)
    {init();};

    /// \brief Constructor from a HDF5 filepath.
    /// @param hdf5_filename Filename of the HDF5 to use for construction.
    /// @param grp HDF5 file group path where the nuSQUIDS object is save.
    /// \details Reads the HDF5 file and construct the associated nuSQUIDS object
    /// restoring all properties as well as the state.
    /// @see ReadStateHDF5
    nuSQUIDS(std::string hdf5_filename, std::string grp = "/") { ReadStateHDF5(hdf5_filename, grp); };

  public:// should this be protected? should this by private? should this exist?
    /************************************************************************************
     * INITIALIZERS
    *************************************************************************************/
    /// \brief Multiple energy mode initializer.
    /// @param Emin Minimum neutrino energy [GeV].
    /// @param Emax Maximum neutirno energy [GeV].
    /// @param Esize Number of energy nodes.
    /// @param numneu Number of neutrino flavors.
    /// @param NT NeutrinoType: neutrino,antineutrino, or both (simultaneous solution).
    /// @param elogscale Sets the energy scale to be logarithmic
    /// @param iinteraction Sets the neutrino noncoherent neutrino interactions on.
    /// \details By default the energy scale is logarithmic and interactions are turn off.
    /// \warning When interactions are present interpolation is performed to precalculate the neutrino
    /// cross section which make take considertable time depending on the energy grid.
    /// @see init
    void Init(double Emin,double Emax,unsigned int Esize,unsigned int numneu_,NeutrinoType NT_ = both,
      bool elogscale_ = true,bool iinteraction_ = false){
      iinteraction = iinteraction_;
      elogscale = elogscale_;
      numneu = numneu_;
      NT = NT_;
      init(Emin,Emax,Esize);
    };

    /// \brief Single energy mode initializer.
    /// @param numneu Number of neutrino flavors.
    /// @param NT NeutrinoType: neutrino or antineutrino.
    /// \warning Interactions are not possible in the single energy mode, nor is simultaneous
    /// neutrino-antineutrino solution (both) possible.
    /// \details Constructors projectors and initializes Hamiltonian.
    void Init(int numneu_, NeutrinoType NT_ = neutrino){
      iinteraction = false;
      elogscale = false;
      numneu = numneu_;
      NT = NT_;
      init();
    };

    /// \brief Initializer from a HDF5 filepath.
    /// @param hdf5_filename Filename of the HDF5 to use for construction.
    /// @param grp HDF5 file group path where the nuSQUIDS object is save.
    /// \details Reads the HDF5 file and construct the associated nuSQUIDS object
    /// restoring all properties as well as the state.
    /// @see ReadStateHDF5
    void Init(std::string hdf5_filename, std::string grp = "/") { ReadStateHDF5(hdf5_filename, grp); };
  protected:
    /************************************************************************************
     * PHYSICS FUNCTIONS - SUPER IMPORTANT
    *************************************************************************************/

    /// \brief Returns the time independent part of the Hamiltonian
    /// @param E Neutrino energy [eV]
    /// @param irho Density matrix equation index.
    /// \details \c irho is the index of the equation on which the Hamiltonian
    /// is used. When not in NeutrinoType == \c both \irho is only 0, but if 
    /// neutrinos and antineutrinos are solved simultaneously, then 0 is for
    /// neutrinos and 1 for antineutrinos.
    SU_vector H0(double E, unsigned int irho) const;

    /// \brief Returns the time dependent part of the Hamiltonian at an energy node \c ie
    /// @param ie Energy node
    /// @param irho Density matrix equation index.
    /// @see H0 for convention on \c irho.
    virtual SU_vector HI(unsigned int ie, unsigned int irho) const;

    /// \brief Returns noncoherent interaction term in the Hamiltonian at node \c ie
    /// @param ie Energy node
    /// @param irho Density matrix equation index.
    /// @see H0 for convention on \c irho.
    virtual SU_vector GammaRho(unsigned int ei, unsigned int irho) const;

    /// \brief Returns the scalar quantity decay rate at an energy node \c ie
    /// @param ie Energy node
    /// @param iscalar scalar equation index.
    /// \details When tau regenereation is considered \c iscalar = 0 corresponds
    /// to the tau flux.
    virtual double GammaScalar(unsigned int ei, unsigned int iscalar) const;

    /// \brief Returns interaction of the density matrix at an energy node \c ie
    /// @param ie Energy node
    /// @param irho Density matrix equation index.
    /// @see H0 for convention on \c irho.
    virtual SU_vector InteractionsRho(unsigned int ei, unsigned int irho) const;

    /// \brief Returns scalar interactions at an energy node \c ie
    /// @param ie Energy node
    /// @param iscalar scalar equation index.
    /// @param irho Density matrix equation index.
    /// @see GammaScalar for convention on \c iscalar.
    virtual double InteractionsScalar(unsigned int ei, unsigned int iscalar) const;
  private:
    /// \brief SQuIDS signature of HI
    SU_vector HI(unsigned int ie, unsigned int irho, double x) const {return HI(ie,irho);};
    /// \brief SQuIDS signature of GammaRho
    SU_vector GammaRho(unsigned int ei, unsigned int irho, double x) const {return GammaRho(ei,irho);};
    /// \brief SQuIDS signature of GammaScalar
    double GammaScalar(unsigned int ei, unsigned int iscalar,double x) const {return GammaScalar(ei,iscalar);};
    /// \brief SQuIDS signature of InteractionsRho
    SU_vector InteractionsRho(unsigned int ei, unsigned int irho, double x) const {return InteractionsRho(ei,irho);};
    /// \brief SQuIDS signature of InteractionsScalar
    double InteractionsScalar(unsigned int ei, unsigned int irho, double x) const {return InteractionsScalar(ei,irho);};
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
    void Set_initial_state(marray<double,1> ini_state, std::string basis = "flavor");

    /// \brief Sets the initial state in the multiple energy mode when doing either neutrino or antineutrino only.
    /// @param ini_state Initial neutrino state.
    /// @param basis Representation of the neutrino state either flavor or mass.
    /// \details \c ini_state first dimension has length equal to \c ne (number of energy nodes), while
    /// the second dimension has length equal to \c numneu (number of flavors). If the basis is
    /// flavor then the entries are interpret as nu_e, nu_mu, nu_tau, nu_sterile_1, ..., nu_sterile_n,
    /// while if the mass basis is used then the first entries correspond to the active
    /// mass eigenstates.
    void Set_initial_state(marray<double,2> ini_state, std::string basis = "flavor");

    /// \brief Sets the initial state in the multiple energy mode when doing either neutrino or antineutrino only.
    /// @param ini_state Initial neutrino state.
    /// @param basis Representation of the neutrino state either flavor or mass.
    /// \details \c ini_state first dimension has length equal to \c ne (number of energy nodes), while
    /// the second dimension has length equal to \c numneu (number of flavors). If the basis is
    /// flavor then the entries are interpret as nu_e, nu_mu, nu_tau, nu_sterile_1, ..., nu_sterile_n,
    /// while if the mass basis is used then the first entries correspond to the active
    /// mass eigenstates.
    void Set_initial_state(marray<double,3> ini_state, std::string basis = "flavor");

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
    double EvalMassAtNode(unsigned int,unsigned int,unsigned int rho = 0) const;
    /// \brief Returns the flavor composition at a given node.
    double EvalFlavorAtNode(unsigned int,unsigned int,unsigned int rho = 0) const;

    /// \brief Returns the mass composition at a given energy in the multiple energy mode.
    double EvalMass(unsigned int,double,unsigned int rho = 0) const;
    /// \brief Returns the flavor composition at a given energy in the multiple energy mode.
    double EvalFlavor(unsigned int,double,unsigned int rho = 0) const;
    /// \brief Returns the mass composition in the single energy mode.
    /// @param flv Neutrino flavor.
    double EvalMass(unsigned int flv) const;
    /// \brief Returns the flavor composition in the single energy mode.
    /// @param flv Neutrino flavor.
    double EvalFlavor(unsigned int flv) const;

    /// \brief Toggles tau regeneration on and off.
    /// @param opt If \c true tau regeneration will be considered.
    void Set_TauRegeneration(bool opt);

    /// \brief Toggles the progress bar printing on and off
    /// @param opt If \c true a progress bar will be printed.
    void Set_ProgressBar(bool opt);

    /// \brief Returns the energy nodes values.
    marray<double,1> GetERange() const;

    /// \brief Returns number of energy nodes.
    size_t GetNumE() const;

    /// \brief Returns the number of neutrino flavors.
    unsigned int GetNumNeu() const;

    /// \brief Return the Hamiltonian.
    SU_vector GetHamiltonian(std::shared_ptr<Track> track, double E, unsigned int rho = 0);
    /// \brief Returns the state
    SU_vector GetState(unsigned int,unsigned int rho = 0) const;
    /// \brief Returns the flavor projector
    SU_vector GetFlavorProj(unsigned int,unsigned int rho = 0) const;
    /// \brief Returns the mass projector
    SU_vector GetMassProj(unsigned int,unsigned int rho = 0) const;

    /// \brief Returns the trajectory object.
    std::shared_ptr<Track> GetTrack() const;
    /// \brief Returns the body object.
    std::shared_ptr<Body> GetBody() const;

    /// \brief Writes the object into an HDF5 file.
    /// @param hdf5_filename Filename of the HDF5 to use for construction.
    /// @param group Path to the group where the nuSQUIDS content will be saved.
    /// @param save_cross_sections If \c true the cross section tables will also be saved.
    /// @param cross_section_grp_loc Path to the group where the cross section will be saved.
    /// \details By default all contents are saved to the \c root of the HDF5 file
    /// if the user wants a different location it can use \c group and \c save_cross_sections to
    /// change the nuSQUIDS location and the cross section information location. Furthermore, the 
    /// \c save_cross_sections flag, which defaults to \c false decides if the cross section
    /// will be saved. Finally, it calls AddToWriteHDF5() to enable user defined parameters
    /// to be saved.
    /// @see AddToWriteHDF5
    /// @see ReadStateHDF5
    void WriteStateHDF5(std::string hdf5_filename,std::string group = "/",bool save_cross_sections = true, std::string cross_section_grp_loc = "") const;
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

    /// \brief Sets mixing parameter.
    /// @param mp Mixing parameter
    /// @param val Value of the parameter
    /// \details Mixing parameters can correspond to mixing angles, CP phases, or 
    /// square mass differences.
    void Set(MixingParameter mp,double val);

    /// \brief Sets the mixing parameters to default.
    void Set_MixingParametersToDefault();

    /// \brief Sets the basis on which the evolution will be perfomed.
    /// @param basis Evolution basis can be either \c mass or \c interaction.
    /// \details By detault the solution is perfomed in the interaction basis, which
    /// is recommended for most problems. One might consider the user of the mass
    /// basis when collective effects, defined by new
    /// interactions introduce phases that cannot be cancelled as one goes to the
    /// interaction basis.
    void Set_Basis(BASIS basis);
};

/**
 * The following class provides functionalities
 * for atmospheric neutrino experiments
 * where a collection of trayectories is explored.
 */

//template<typename BaseType = nuSQUIDS, typename = typename std::enable_if<std::is_base_of<nuSQUIDS,BaseType>>::type >
class nuSQUIDSAtm {
  private:
    bool progressbar = false;
    double LinInter(double,double,double,double,double) const;
  protected:
    bool iinistate;
    bool inusquidsatm;
    marray<double,1> costh_array;
    marray<double,1> enu_array;
    marray<double,1> log_enu_array;
    std::vector<nuSQUIDS> nusq_array;

    std::shared_ptr<EarthAtm> earth_atm;
    std::vector<std::shared_ptr<EarthAtm::Track>> track_array;
  public:
    nuSQUIDSAtm(double,double,int,double,double,int,
                int numneu,NeutrinoType NT = both,
                bool elogscale = true, bool iinteraction = false);
    nuSQUIDSAtm(std::string str) {ReadStateHDF5(str);};

    void Set_initial_state(marray<double,3>, std::string basis = "flavor");
    void Set_initial_state(marray<double,4>, std::string basis = "flavor");

    void EvolveState();
    void Set_TauRegeneration(bool);

    double EvalFlavor(unsigned int,double,double,unsigned int rho = 0) const;

    void WriteStateHDF5(std::string) const;
    void ReadStateHDF5(std::string);
    void Set_MixingParametersToDefault(void);
    void Set(MixingParameter,double);
    void Set_abs_error(double);
    void Set_rel_error(double);
    const Const units;

    void Set_ProgressBar(bool);

    size_t GetNumE() const;
    size_t GetNumCos() const;
    marray<double,1> GetERange() const;
    marray<double,1> GetCosthRange() const;
};


} // close namespace
#endif
