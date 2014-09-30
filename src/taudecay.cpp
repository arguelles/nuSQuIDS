#include "taudecay.h"

#define SQR(x)      ((x)*(x))                        // x^2

// decay formulaes
namespace nusquids{

double TauDecaySpectra::TauDecayToLepton(double E_tau,double E_nu){
    double z = E_nu/E_tau;
    double g0 = 5.0/3.0-3.0*pow(z,2.0)+(4.0/3.0)*pow(z,3.0);
    double g1 = 1.0/3.0-3.0*pow(z,2.0)+(8.0/3.0)*pow(z,3.0);

    return g0+TauPolarization*g1;
}

double TauDecaySpectra::TauDecayToPion(double E_tau,double E_nu){
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

double TauDecaySpectra::TauDecayToRho(double E_tau,double E_nu){
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

double TauDecaySpectra::TauDecayToA1(double E_tau,double E_nu){
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

double TauDecaySpectra::TauDecayToHadron(double E_tau,double E_nu){
    double z = E_nu/E_tau;

    double g0 = 0.0;
    if (0.3 - z > 0.0) {
        g0 = 1.0/0.3;
    }

    double g1 = 0.0;

    return g0+TauPolarization*g1;
}

double TauDecaySpectra::TauDecayToAll(double E_tau, double E_nu){
    double decay_spectra = 0.0;

    decay_spectra += 2.0*BrLepton*TauDecayToLepton(E_tau,E_nu);
    decay_spectra += BrPion*TauDecayToPion(E_tau,E_nu);
    decay_spectra += BrRho*TauDecayToRho(E_tau,E_nu);
    decay_spectra += BrRA1*TauDecayToA1(E_tau,E_nu);
    decay_spectra += BrHadron*TauDecayToHadron(E_tau,E_nu);

    return decay_spectra;
}

// define tau decay object

TauDecaySpectra::TauDecaySpectra(){};

void TauDecaySpectra::SetParameters(void){
  TauPolarization = -1.0;
  RPion = SQR(0.07856),RRho = SQR(0.43335),RA1 = SQR(0.70913);
  BrLepton = 0.18,BrPion = 0.12,BrRho = 0.26,BrRA1 = 0.13,BrHadron = 0.13;
}

TauDecaySpectra::TauDecaySpectra(double Emin_in,double Emax_in,int div_in){
  Init(Emin_in,Emax_in,div_in);
};

void TauDecaySpectra::Init(double Emin_in,double Emax_in,int div_in){
        SetParameters();
        Emin = Emin_in;
        Emax = Emax_in;
        div = div_in;

        std::vector<double> E_range_GeV = logspace(Emin/1.0e9,Emax/1.0e9,div);
        int e_size = E_range_GeV.size();

        for (int e1 = 0 ; e1 < e_size ; e1 ++){
            double Enu1 = E_range_GeV[e1];
            for (int e2 = 0 ; e2 < e_size ; e2 ++){
                double Enu2 = E_range_GeV[e2];
                // create rows
                Row dNdEle_All;
                Row dNdEle_Lep;
                Row dNdEnu_All;
                Row dNdEnu_Lep;

                // save energies in table
                dNdEnu_All.push_back(Enu1);
                dNdEnu_Lep.push_back(Enu1);
                dNdEnu_All.push_back(Enu2);
                dNdEnu_Lep.push_back(Enu2);

                dNdEle_All.push_back(Enu1);
                dNdEle_Lep.push_back(Enu1);
                dNdEle_All.push_back(Enu2);
                dNdEle_Lep.push_back(Enu2);

                // save spectra
                dNdEle_All.push_back(TauDecayToAll(Enu1,Enu2)*Enu2/(Enu1*Enu1));
                dNdEle_Lep.push_back(BrLepton*TauDecayToLepton(Enu1,Enu2)*Enu2/(Enu1*Enu1));

                dNdEnu_All.push_back(TauDecayToAll(Enu1,Enu2)/Enu1);
                dNdEnu_Lep.push_back(BrLepton*TauDecayToLepton(Enu1,Enu2)/Enu1);


                // append row
                dNdEnu_All_tbl.push_back(dNdEnu_All);
                dNdEnu_Lep_tbl.push_back(dNdEnu_Lep);

                dNdEle_All_tbl.push_back(dNdEle_All);
                dNdEle_Lep_tbl.push_back(dNdEle_Lep);
            }
        }

};

// tau decay spectra returned in units of [GeV^-1]

double TauDecaySpectra::dNdEnu_All(int i_enu, int i_ele){
    int ii = i_enu*(div+1) + (i_ele);
    return dNdEnu_All_tbl[ii][2];
};

double TauDecaySpectra::dNdEnu_Lep(int i_enu, int i_ele){
    int ii = i_enu*(div+1) + i_ele;
    return (double) dNdEnu_Lep_tbl[ii][2];
};

double TauDecaySpectra::dNdEle_All(int i_enu, int i_ele){
        int ii = i_enu*(div+1) + (i_ele);
        return dNdEle_All_tbl[ii][2];
};

double TauDecaySpectra::dNdEle_Lep(int i_enu, int i_ele){
        int ii = i_enu*(div+1) + i_ele;
        return (double) dNdEle_Lep_tbl[ii][2];
};

}// close namespace
