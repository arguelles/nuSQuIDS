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


#include <nuSQuIDS/taudecay.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#define SQR(x)      ((x)*(x))                        // x^2

// decay formulaes
namespace nusquids{

// in the following functions E_tau refers to the
// decaying tau energy, while E_nu refers to the secondary tau

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

void TauDecaySpectra::SetParameters(bool neutrino_type){
  TauPolarization = -1.0; // due to CP symmetry both tau and antitau spectra are the same
  RPion = SQR(0.07856); RRho = SQR(0.43335); RA1 = SQR(0.70913);
  BrLepton = 0.18; BrPion = 0.12; BrRho = 0.26; BrRA1 = 0.13; BrHadron = 0.13;
}

TauDecaySpectra::TauDecaySpectra(marray<double,1> E_range):TauDecaySpectra(){
  Init(E_range);
}

void TauDecaySpectra::Init(marray<double,1> E_range){
  SetParameters(true);
  double GeV = 1.0e9;
  unsigned int e_size = E_range.size();
  unsigned int neutrino_number = 2;

  dNdEnu_All_tbl.resize(std::vector<size_t>{neutrino_number,e_size,e_size});
  dNdEnu_Lep_tbl.resize(std::vector<size_t>{neutrino_number,e_size,e_size});
  dNdEle_All_tbl.resize(std::vector<size_t>{neutrino_number,e_size,e_size});
  dNdEle_Lep_tbl.resize(std::vector<size_t>{neutrino_number,e_size,e_size});

  for(unsigned int neutype = 0; neutype < neutrino_number; neutype++){
    TauPolarization = -1.0; // due to CP symmetry both tau and antitau spectra are the same
    for(unsigned int e1 = 0 ; e1 < e_size ; e1 ++){
        double Enu1 = E_range[e1]/GeV; // tau energy
        for(unsigned int e2 = 0 ; e2 < e_size ; e2 ++){
            double Enu2 = E_range[e2]/GeV; // tau neutrino energy
            // save spectra
            dNdEle_All_tbl[neutype][e1][e2] = TauDecayToAll(Enu1,Enu2)*Enu2/(Enu1*Enu1);
            dNdEle_Lep_tbl[neutype][e1][e2] = BrLepton*TauDecayToLepton(Enu1,Enu2)*Enu2/(Enu1*Enu1);

            dNdEnu_All_tbl[neutype][e1][e2] = TauDecayToAll(Enu1,Enu2)/Enu1;
            dNdEnu_Lep_tbl[neutype][e1][e2] = BrLepton*TauDecayToLepton(Enu1,Enu2)/Enu1;
        }
    }
  }
}

// tau decay spectra returned in units of [GeV^-1]

double TauDecaySpectra::dNdEnu_All(unsigned int i_enu,unsigned int i_ele, unsigned int neutrino_type) const{
  if(i_enu >= dNdEnu_All_tbl.extent(1) or i_ele >= dNdEnu_All_tbl.extent(2))
    throw std::runtime_error("TauDecaySpectra:dNdEnu_All in valid index choices");
  return dNdEnu_All_tbl[neutrino_type][i_enu][i_ele];
}

double TauDecaySpectra::dNdEnu_Lep(unsigned int i_enu,unsigned int i_ele, unsigned int neutrino_type) const{
  if(i_enu >= dNdEnu_Lep_tbl.extent(1) or i_ele >= dNdEnu_Lep_tbl.extent(2))
    throw std::runtime_error("TauDecaySpectra:dNdEnu_Lep in valid index choices");
  return dNdEnu_Lep_tbl[neutrino_type][i_enu][i_ele];
}

double TauDecaySpectra::dNdEle_All(unsigned int i_enu,unsigned int i_ele, unsigned int neutrino_type) const{
  if(i_enu >= dNdEle_All_tbl.extent(1) or i_ele >= dNdEle_All_tbl.extent(2))
    throw std::runtime_error("TauDecaySpectra:dNdEle_All in valid index choices");
  return dNdEle_All_tbl[neutrino_type][i_enu][i_ele];
}

double TauDecaySpectra::dNdEle_Lep(unsigned int i_enu,unsigned int i_ele, unsigned int neutrino_type) const{
  if(i_enu >= dNdEle_Lep_tbl.extent(1) or i_ele >= dNdEle_Lep_tbl.extent(2))
    throw std::runtime_error("TauDecaySpectra:dNdEle_Lepin valid index choices");
  return dNdEle_Lep_tbl[neutrino_type][i_enu][i_ele];
}

}// close namespace
