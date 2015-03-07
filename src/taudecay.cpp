#include "taudecay.h"

#define SQR(x)      ((x)*(x))                        // x^2

// decay formulaes
namespace nusquids{

double TauDecaySpectra::TauDecayToLepton(double E_tau,double E_nu) const{
    double z = E_nu/E_tau;
    double g0 = 5.0/3.0-3.0*pow(z,2.0)+(4.0/3.0)*pow(z,3.0);
    double g1 = 1.0/3.0-3.0*pow(z,2.0)+(8.0/3.0)*pow(z,3.0);

    return g0+TauPolarization*g1;
}

double TauDecaySpectra::TauDecayToPion(double E_tau,double E_nu) const{
    double z = E_nu/E_tau;

    double g0 = 0.0;
    if (1.0 - RPion - z > 0.0) {
        g0 = 1.0/(1.0 - RPion);
    }

    double g1 = 0.0;
    if (1.0 - RPion - z > 0.0) {
        g1 = -(2.0*z-1.0+RPion)/pow(1.0-RPion,2.0);
    }

    return g0+TauPolarization*g1;
}

double TauDecaySpectra::TauDecayToRho(double E_tau,double E_nu) const{
    double z = E_nu/E_tau;

    double g0 = 0.0;
    if (1.0 - RRho - z > 0.0) {
        g0 = 1.0/(1.0 - RRho);
    }

    double g1 = 0.0;
    if (1.0 - RRho - z > 0.0) {
        g1 = -((2.0*z-1.0+RRho)/(1.0-RRho))*((1.0-2.0*RRho)/(1.0+2.0*RRho));
    }

    return g0+TauPolarization*g1;
}

double TauDecaySpectra::TauDecayToA1(double E_tau,double E_nu) const{
    double z = E_nu/E_tau;

    double g0 = 0.0;
    if (1.0 - RA1 - z > 0.0) {
        g0 = 1.0/(1.0 - RA1);
    }

    double g1 = 0.0;
    if (1.0 - RA1 - z > 0.0) {
        g1 = -((2.0*z-1.0+RA1)/(1.0-RA1))*((1.0-2.0*RA1)/(1.0+2.0*RA1));
    }

    return g0+TauPolarization*g1;
}

double TauDecaySpectra::TauDecayToHadron(double E_tau,double E_nu) const{
    double z = E_nu/E_tau;

    double g0 = 0.0;
    if (0.3 - z > 0.0) {
        g0 = 1.0/0.3;
    }

    double g1 = 0.0;

    return g0+TauPolarization*g1;
}

double TauDecaySpectra::TauDecayToAll(double E_tau, double E_nu) const{
    double decay_spectra = 0.0;

    decay_spectra += 2.0*BrLepton*TauDecayToLepton(E_tau,E_nu);
    decay_spectra += BrPion*TauDecayToPion(E_tau,E_nu);
    decay_spectra += BrRho*TauDecayToRho(E_tau,E_nu);
    decay_spectra += BrRA1*TauDecayToA1(E_tau,E_nu);
    decay_spectra += BrHadron*TauDecayToHadron(E_tau,E_nu);

    return decay_spectra;
}

// define tau decay object

TauDecaySpectra::TauDecaySpectra(){}

void TauDecaySpectra::SetParameters(){
  TauPolarization = -1.0;
  RPion = SQR(0.07856),RRho = SQR(0.43335),RA1 = SQR(0.70913);
  BrLepton = 0.18,BrPion = 0.12,BrRho = 0.26,BrRA1 = 0.13,BrHadron = 0.13;
}

TauDecaySpectra::TauDecaySpectra(double Emin_in,double Emax_in,unsigned int div_in){
  Init(Emin_in,Emax_in,div_in);
}

void TauDecaySpectra::Init(double Emin_in,double Emax_in,unsigned int div_in){
        SetParameters();
        Emin = Emin_in;
        Emax = Emax_in;
        div = div_in;

        marray<double,1> E_range_GeV = logspace(Emin/1.0e9,Emax/1.0e9,div);
        unsigned int e_size = E_range_GeV.size();

        dNdEnu_All_tbl.resize(std::vector<size_t>{e_size,e_size});
        dNdEnu_Lep_tbl.resize(std::vector<size_t>{e_size,e_size});
        dNdEle_All_tbl.resize(std::vector<size_t>{e_size,e_size});
        dNdEle_Lep_tbl.resize(std::vector<size_t>{e_size,e_size});

        for (unsigned int e1 = 0 ; e1 < e_size ; e1 ++){
            double Enu1 = E_range_GeV[e1];
            for (unsigned int e2 = 0 ; e2 < e_size ; e2 ++){
                double Enu2 = E_range_GeV[e2];
                // save spectra
                dNdEle_All_tbl[e1][e2] = TauDecayToAll(Enu1,Enu2)*Enu2/(Enu1*Enu1);
                dNdEle_Lep_tbl[e1][e2] = BrLepton*TauDecayToLepton(Enu1,Enu2)*Enu2/(Enu1*Enu1);

                dNdEnu_All_tbl[e1][e2] = TauDecayToAll(Enu1,Enu2)/Enu1;
                dNdEnu_Lep_tbl[e1][e2] = BrLepton*TauDecayToLepton(Enu1,Enu2)/Enu1;
            }
        }

}

// tau decay spectra returned in units of [GeV^-1]

double TauDecaySpectra::dNdEnu_All(unsigned int i_enu,unsigned int i_ele) const{
    return dNdEnu_All_tbl[i_enu][i_ele];
}

double TauDecaySpectra::dNdEnu_Lep(unsigned int i_enu,unsigned int i_ele) const{
    return dNdEnu_Lep_tbl[i_enu][i_ele];
}

double TauDecaySpectra::dNdEle_All(unsigned int i_enu,unsigned int i_ele) const{
    return dNdEle_All_tbl[i_enu][i_ele];
}

double TauDecaySpectra::dNdEle_Lep(unsigned int i_enu,unsigned int i_ele) const{
    return dNdEle_Lep_tbl[i_enu][i_ele];
}

}// close namespace
