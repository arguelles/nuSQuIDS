#ifndef __nuSQUID_H
#define __nuSQUID_H

#include "body.h"
#include "xsections.h"
#include "taudecay.h"

#include <SQUIDS/SQUIDS.h>
#include <memory>
#include <map>
#include <stdexcept>

//#include "H5Exception.h"
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

typedef std::vector<double> array1D;
typedef std::vector<std::vector<double> > array2D;
typedef std::vector<std::vector<std::vector<double> > > array3D;
typedef std::vector<std::vector<std::vector<std::vector<double> > > > array4D;

enum MixingParameter {
  TH12 = 0,
  TH13 = 1, TH23 = 2,
  TH14 = 3, TH24 = 4, TH34 = 5,
  TH15 = 6, TH25 = 7, TH35 = 8, TH45 = 9,
  TH16 = 10, TH26 = 11, TH36 = 12, TH46 = 13, TH56 = 14,
  DELTA1 = 15, DELTA2 = 16, DELTA3 = 17,
  DM21SQ = 18, DM31SQ = 19, DM41SQ = 20, DM51SQ = 21, DM61SQ = 22
                    };

enum mode { neutrino, antineutrino, both };

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
  protected:
    BASIS basis = interaction;
    int numneu;
    int ne = nx;

    double GetNucleonNumber();
    void UpdateInteractions();

    array1D E_range;
    array1D delE;

    NeutrinoCrossSections ncs;
    array4D dNdE_CC,dNdE_NC;
    array3D invlen_NC,invlen_CC,invlen_INT;
    array3D sigma_CC,sigma_NC;

    TauDecaySpectra tdc;
    array1D invlen_tau;
    array2D dNdE_tau_all,dNdE_tau_lep;
    double taubr_lep,tau_lifetime,tau_mass;
    double tau_reg_scale;

    std::shared_ptr<Body> body;
    std::shared_ptr<Track> track;

    SU_vector DM2;
    std::vector<SU_vector> H0_array;

    std::vector<SU_vector> b0_proj;
    std::vector<std::vector<SU_vector> > b1_proj;
    std::vector<std::vector<std::vector<SU_vector> > > evol_b0_proj;
    std::vector<std::vector<std::vector<SU_vector> > > evol_b1_proj;

    std::vector<std::vector<SU_vector> > potential_array;

    void EvolveProjectors(double t);
    void ConvertTauIntoNuTau(void);

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
    std::string NT = "both"; // {"neutrino","antineutrino","both"}

    void iniProjectors();
    void SetIniFlavorProyectors();
    void iniH0();
    void AntineutrinoCPFix(int);
    void SetBodyTrack(int,int,double*,int,double*);

     void init(double,double,int,bool initialize_intereractions = true);
     void init();
     void InitializeInteractionVectors(void);
     void InitializeInteractions(void);
     void Set_Initial_Time(void);
     void SetScalarsToZero(void);
     void ProgressBar(void);
  public:
     virtual void PreDerive(double);
     virtual void AddToPreDerive(double){};
     const Const units;
     // initializers
     nuSQUIDS(){};
     nuSQUIDS(double Emin,double Emax,int Esize,int numneu,std::string NT = "both",
         bool elogscale = true,bool iinteraction = false):
     iinteraction(iinteraction),elogscale(elogscale),numneu(numneu),NT(NT)
     {init(Emin,Emax,Esize);};

     void Init(double Emin,double Emax,int Esize,int numneu_,std::string NT_ = "both",
         bool elogscale_ = true,bool iinteraction_ = false){
       iinteraction = iinteraction_;
       elogscale = elogscale_;
       numneu = numneu_;
       NT = NT_;
       init(Emin,Emax,Esize);
     };

     nuSQUIDS(int numneu, std::string NT = "neutrino"):
     iinteraction(false),elogscale(false),numneu(numneu),NT(NT)
     {init();};

     void Init(int numneu_, std::string NT_ = "neutrino"){
      iinteraction = false;
      elogscale = false;
      numneu = numneu_;
      NT = NT_;
      init();
     };

     nuSQUIDS(std::string in_hdf5, std::string grp = "/") { ReadStateHDF5(in_hdf5, grp); };
     void Init(std::string in_hdf5, std::string grp = "/") { ReadStateHDF5(in_hdf5, grp); };
     // physics functions
     SU_vector H0(double);

     virtual SU_vector HI(int);
     virtual SU_vector HI(int ei,double x){return tunit*HI(ei);};

     virtual SU_vector GammaRho(int);
     virtual SU_vector GammaRho(int ei,double x){return tunit*GammaRho(ei);};
     virtual double GammaScalar(int);
     virtual double GammaScalar(int ei,double x){return tunit*GammaScalar(ei);};

     virtual SU_vector InteractionsRho(int);
     virtual SU_vector InteractionsRho(int ei,double x){return tunit*InteractionsRho(ei);};
     virtual double InteractionsScalar(int);
     virtual double InteractionsScalar(int ei,double x){return tunit*InteractionsScalar(ei);};
     // interface
     void Set_initial_state(array1D, std::string basis = "flavor");
     void Set_initial_state(array2D, std::string basis = "flavor");
     void Set_initial_state(array3D, std::string basis = "flavor");

     void Set_Body(std::shared_ptr<Body>);
     void Set_Track(std::shared_ptr<Track>);

     void Set_E(double);

     void EvolveState(void);

     double EvalMassAtNode(int,int,int rho = 0);
     double EvalFlavorAtNode(int,int,int rho = 0);

     double EvalMass(int,double,int rho = 0);
     double EvalFlavor(int,double,int rho = 0);
     double EvalMass(int);
     double EvalFlavor(int);

     void Set_TauRegeneration(bool);
     void Set_ProgressBar(bool);

     array1D GetERange(void);
     size_t GetNumE(void);
     int GetNumNeu(void);
     SU_vector GetHamiltonian(std::shared_ptr<Track> track, double E, int rho = 0);
     SU_vector GetState(int,int rho = 0);
     SU_vector GetFlavorProj(int,int rho = 0);
     SU_vector GetMassProj(int,int rho = 0);

     std::shared_ptr<Track> GetTrack(void);
     std::shared_ptr<Body> GetBody(void);

     void WriteStateHDF5(std::string,std::string group = "/",bool save_cross_sections = true, std::string cross_section_grp_loc = "");
     void ReadStateHDF5(std::string,std::string group = "/", std::string cross_section_grp_loc = "");

     void Set(MixingParameter,double);
     void Set_MixingParametersToDefault(void);

     void Set_Basis(BASIS);

    // virtual ~nuSQUIDS(void);
};

/**
 * The following class provides functionalities
 * for atmospheric neutrino experiments
 * where a collection of trayectories is explored.
 */

class nuSQUIDSAtm {
  private:
    bool progressbar = false;
    double LinInter(double,double,double,double,double) const;
  protected:
    bool iinistate;
    bool inusquidsatm;
    std::vector<double> costh_array;
    std::vector<double> enu_array;
    std::vector<double> log_enu_array;
    std::vector<nuSQUIDS> nusq_array;

    std::shared_ptr<EarthAtm> earth_atm;
    std::vector<std::shared_ptr<EarthAtm::Track>> track_array;
  public:
    nuSQUIDSAtm(double,double,int,double,double,int,
                int numneu,std::string NT = "both",
                bool elogscale = true, bool iinteraction = false);
    nuSQUIDSAtm(std::string str) {ReadStateHDF5(str);};

    void Set_initial_state(array3D, std::string basis = "flavor");
    void Set_initial_state(array4D, std::string basis = "flavor");

    void EvolveState(void);
    void Set_TauRegeneration(bool);

    double EvalFlavor(int,double,double,int rho = 0);

    void WriteStateHDF5(std::string);
    void ReadStateHDF5(std::string);
    void Set_MixingParametersToDefault(void);
    void Set(MixingParameter,double);
    void Set_abs_error(double);
    void Set_rel_error(double);
    const Const units;

    void Set_ProgressBar(bool);

    size_t GetNumE(void);
    size_t GetNumCos(void);
    array1D GetERange(void);
    array1D GetCosthRange(void);
};


} // close namespace
#endif
