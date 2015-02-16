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

class nuSQUIDS: public SQUIDS {
  public:
    enum NeutrinoType {
      neutrino=0b01,
      antineutrino=0b10,
      both=0b11
    };
  protected:
    BASIS basis = interaction;
    unsigned int numneu;
    unsigned int ne = nx;

    double GetNucleonNumber();
    void UpdateInteractions();

    marray<double,1> E_range;
    marray<double,1> delE;

    //void Set_Initial_Time();

    NeutrinoCrossSections ncs;
    marray<double,4> dNdE_CC,dNdE_NC;
    marray<double,3> invlen_NC,invlen_CC,invlen_INT;
    marray<double,3> sigma_CC,sigma_NC;

    TauDecaySpectra tdc;
    marray<double,1> invlen_tau;
    marray<double,2> dNdE_tau_all,dNdE_tau_lep;
    double taubr_lep,tau_lifetime,tau_mass;
    double tau_reg_scale;

    std::shared_ptr<Body> body;
    std::shared_ptr<Track> track;

    SU_vector DM2;
    marray<SU_vector,1> H0_array;

    marray<SU_vector,1> b0_proj;
    marray<SU_vector,2> b1_proj;
    marray<SU_vector,3> evol_b0_proj;
    marray<SU_vector,3> evol_b1_proj;

    marray<SU_vector,2> potential_array;

    void EvolveProjectors(double t);
    void ConvertTauIntoNuTau();

    // bool requirements
    bool inusquids = false;
    bool ibody = false;
    bool ienergy = false;
    bool itrack = false;
    bool ioscpar = false;
    bool istate = false;
    bool iinteraction = false;
    bool elogscale = true;
    bool tauregeneration = false;
    bool progressbar = false;
    int progressbar_count = 0;
    int progressbar_loop = 100;
    // neutrino type
    NeutrinoType NT = both;

    void iniProjectors();
    void SetIniFlavorProyectors();
    void iniH0();
    void AntineutrinoCPFix(unsigned int irho);
    void SetBodyTrack(int,int,double*,int,double*);

     void init(double,double,unsigned int,bool initialize_intereractions = true, double xini = 0.0);
     void init(double xini = 0.0);
     void InitializeInteractionVectors();
     void InitializeInteractions();
     void SetScalarsToZero();
     void ProgressBar() const;
  public:
     virtual void PreDerive(double);
     virtual void AddToPreDerive(double){};
     const Const units;
     // initializers
     nuSQUIDS(){};
     nuSQUIDS(double Emin,double Emax,unsigned int Esize,unsigned int numneu,NeutrinoType NT = both,
         bool elogscale = true,bool iinteraction = false):
     iinteraction(iinteraction),elogscale(elogscale),numneu(numneu),NT(NT)
     {init(Emin,Emax,Esize);};

     void Init(double Emin,double Emax,unsigned int Esize,unsigned int numneu_,NeutrinoType NT_ = both,
         bool elogscale_ = true,bool iinteraction_ = false){
       iinteraction = iinteraction_;
       elogscale = elogscale_;
       numneu = numneu_;
       NT = NT_;
       init(Emin,Emax,Esize);
     };

     nuSQUIDS(int numneu, NeutrinoType NT = neutrino):
     iinteraction(false),elogscale(false),numneu(numneu),NT(NT)
     {init();};

     void Init(int numneu_, NeutrinoType NT_ = neutrino){
      iinteraction = false;
      elogscale = false;
      numneu = numneu_;
      NT = NT_;
      init();
     };

     nuSQUIDS(std::string in_hdf5, std::string grp = "/") { ReadStateHDF5(in_hdf5, grp); };
     void Init(std::string in_hdf5, std::string grp = "/") { ReadStateHDF5(in_hdf5, grp); };
     // physics functions
     SU_vector H0(double E, unsigned int irho) const;

     virtual SU_vector HI(unsigned int ie, unsigned int irho) const;
     virtual SU_vector HI(unsigned int ie, unsigned int irho, double x) const {return HI(ie,irho);};

     virtual SU_vector GammaRho(unsigned int ei, unsigned int irho) const;
     virtual SU_vector GammaRho(unsigned int ei, unsigned int irho, double x) const {return GammaRho(ei,irho);};
     virtual double GammaScalar(unsigned int ei, unsigned int iscalar) const;
     virtual double GammaScalar(unsigned int ei, unsigned int iscalar,double x) const {return GammaScalar(ei,iscalar);};

     virtual SU_vector InteractionsRho(unsigned int ei, unsigned int irho) const;
     virtual SU_vector InteractionsRho(unsigned int ei, unsigned int irho, double x) const {return InteractionsRho(ei,irho);};
     virtual double InteractionsScalar(unsigned int ei, unsigned int irho) const;
     virtual double InteractionsScalar(unsigned int ei, unsigned int irho, double x) const {return InteractionsScalar(ei,irho);};
     // interface
     void Set_initial_state(marray<double,1>, std::string basis = "flavor");
     void Set_initial_state(marray<double,2>, std::string basis = "flavor");
     void Set_initial_state(marray<double,3>, std::string basis = "flavor");

     void Set_Body(std::shared_ptr<Body>);
     void Set_Track(std::shared_ptr<Track>);

     void Set_E(double);

     void EvolveState();

     double EvalMassAtNode(unsigned int,unsigned int,unsigned int rho = 0) const;
     double EvalFlavorAtNode(unsigned int,unsigned int,unsigned int rho = 0) const;

     double EvalMass(unsigned int,double,unsigned int rho = 0) const;
     double EvalFlavor(unsigned int,double,unsigned int rho = 0) const;
     double EvalMass(unsigned int) const;
     double EvalFlavor(unsigned int) const;

     void Set_TauRegeneration(bool);
     void Set_ProgressBar(bool);

     marray<double,1> GetERange() const;
     size_t GetNumE() const;
     int GetNumNeu() const;
     SU_vector GetHamiltonian(std::shared_ptr<Track> track, double E, unsigned int rho = 0);
     SU_vector GetState(unsigned int,unsigned int rho = 0) const;
     SU_vector GetFlavorProj(unsigned int,unsigned int rho = 0) const;
     SU_vector GetMassProj(unsigned int,unsigned int rho = 0) const;

     std::shared_ptr<Track> GetTrack() const;
     std::shared_ptr<Body> GetBody() const;

     void WriteStateHDF5(std::string,std::string group = "/",bool save_cross_sections = true, std::string cross_section_grp_loc = "") const;
     virtual void AddToWriteHDF5(hid_t hdf5_loc_id) const;
     void ReadStateHDF5(std::string,std::string group = "/", std::string cross_section_grp_loc = "");
     virtual void AddToReadHDF5(hid_t hdf5_loc_id);

     void Set(MixingParameter,double);
     void Set_MixingParametersToDefault();

     void Set_Basis(BASIS);
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
                int numneu,nuSQUIDS::NeutrinoType NT = nuSQUIDS::both,
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
