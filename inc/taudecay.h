#ifndef __TAUDECAY_H
#define __TAUDECAY_H

#include "tools.h"
#include <string>
#include <cmath>
#include <math.h>

/*
 * FORMULAES implemented out of
 *
 * Tau neutrinos underground
 * PHYSICAL REVIEW D 62 123001
 * SHARADA IYER DUTTA, MARY HALL RENO, AND INA SARCEVIC
 *
 */

// Tau decay spectra

class TauDecaySpectra{
    private:
      double TauPolarization;
      double RPion,RRho,RA1;
      double BrLepton,BrPion,BrRho,BrRA1,BrHadron;

      double TauDecayToLepton(double,double);
      double TauDecayToPion(double,double);
      double TauDecayToRho(double,double);
      double TauDecayToA1(double,double);
      double TauDecayToHadron(double,double);
      double TauDecayToAll(double,double);

      void SetParameters(void);

      Table dNdEnu_All_tbl;
      Table dNdEnu_Lep_tbl;

      Table dNdEle_All_tbl;
      Table dNdEle_Lep_tbl;

      double Emin;
      double Emax;
      int div;
    public :
      TauDecaySpectra();
      TauDecaySpectra(double,double,int);
      void Init(double,double,int);

      double dNdEnu_All(int,int);
      double dNdEnu_Lep(int,int);

      double dNdEle_All(int,int);
      double dNdEle_Lep(int,int);
};

#endif
