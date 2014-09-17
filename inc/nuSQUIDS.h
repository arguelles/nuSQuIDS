#ifndef __nuSQUID_H
#define __nuSQUID_H

#include <body.h>
#include <xsections.h>
#include <taudecay.h>
#include <SQUIDS.h>
#include <memory>
#include <map>

#include "hdf5.h"
#include "hdf5_hl.h"

#define ManualTauReinjection
#define FixCrossSections

//#define UpdateInteractions_DEBUG

typedef vector<double> array1D;
typedef vector<vector<double> > array2D;
typedef vector<vector<vector<double> > > array3D;
typedef vector<vector<vector<vector<double> > > > array4D;

std::map<int,string> param_label_map {
  {0 , "th12"},
  {1 , "th13"}, {2 , "th23"},
  {3 , "th14"}, {4 , "th24"}, {5 , "th34"},
  {6 , "th15"}, {7 , "th25"}, {8 , "th35"}, {9 , "th45"},
  {10 , "th16"}, {11 , "th26"}, {12 , "th36"}, {13 , "th46"}, {14 , "th56"},
  {15 , "delta1"} , {16, "delta2"} , {17, "delta3"},
  {18, "dm21sq"} , {19 , "dm31sq"}, {20,"dm41sq"} ,{ 21 , "dm51sq"} , {22, "dm61sq"}
                                     };

std::map<int,vector<int>> param_label_index {
  {0 , {1,2}},
  {1 , {1,3}}, {2 , {2,3}},
  {3 , {1,4}}, {4 , {2,4}}, {5 , {3,4}},
  {6 , {1,5}}, {7 , {2,5}}, {8 , {3,5}}, {9 , {4,5}},
  {10 , {1,6}},{11 , {2,6}}, {12 , {3,6}}, {13 , {4,6}}, {14 , {5,6}},
  {15 , {1,0}} , {16, {2,0}} , {17, {3,0}},
  {18, {1,0}} , {19 , {2,0}}, {20,{3,0}} ,{ 21 , {4,0}} , {22, {5,0}}
                                       };

class nuSQUIDS: public SQUIDS {
  protected:
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

    vector<SU_vector> b0_proj;
    vector<vector<SU_vector> > b1_proj;
    vector<vector<vector<SU_vector> > > evol_b0_proj;
    vector<vector<vector<SU_vector> > > evol_b1_proj;

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
    // neutrino type
    string NT = "both"; // {"neutrino","antineutrino","both"}

    void iniProyectors();
    void SetIniFlavorProyectors();
    void iniH0();

    void SetBodyTrack(int,int,double*,int,double*);

     void init(double,double,int);
     void init();

     void InitializeInteractions(void);
     void Set_Initial_Time(void);
     void SetScalarsToZero(void);
  public:
     void PreDerive(double);
     const Const units;
     // initializers
     nuSQUIDS(){};
     nuSQUIDS(double Emin,double Emax,int Esize,int numneu,string NT = "both",
         bool elogscale = true,bool iinteraction = false):
     iinteraction(iinteraction),elogscale(elogscale),numneu(numneu),NT(NT)
     {init(Emin,Emax,Esize);};

     void Init(double Emin,double Emax,int Esize,int numneu_,string NT_ = "both",
         bool elogscale_ = true,bool iinteraction_ = false){
       iinteraction = iinteraction_;
       elogscale = elogscale_;
       numneu = numneu_;
       NT = NT_;
       init(Emin,Emax,Esize);
     };

     nuSQUIDS(int numneu, string NT = "neutrino"):
     iinteraction(false),elogscale(false),numneu(numneu),NT(NT)
     {init();};

     void Init(int numneu_, string NT_ = "neutrino"){
      iinteraction = false;
      elogscale = false;
      numneu = numneu_;
      NT = NT_;
      init();
     };

     nuSQUIDS(string in_hdf5) { ReadStateHDF5(in_hdf5); };
     void Init(string in_hdf5) { ReadStateHDF5(in_hdf5); };
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
     void Set_initial_state(array1D, string basis = "flavor");
     void Set_initial_state(array2D, string basis = "flavor");
     void Set_initial_state(array3D, string basis = "flavor");

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

     void Set_nuSQUIDS(string,bool);

     array1D GetERange(void);
     size_t GetNumE(void);
     int GetNumNeu(void);
     SU_vector GetHamiltonian(std::shared_ptr<Track> track, double E, int rho = 0);

     void WriteState(string);
     void ReadState(string);

     std::shared_ptr<Track> GetTrack(void);
     std::shared_ptr<Body> GetBody(void);

     void WriteStateHDF5(string);
     void ReadStateHDF5(string);

    // virtual ~nuSQUIDS(void);
};

#endif
