#ifndef EXOBJ_H
#define EXOBJ_H

#include <unistd.h>
#include <iostream>
#include <fstream>

#include <float.h>
#include <math.h>
#include <complex>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>

namespace nusquids {

// LinearCrossSections
// Simple toy cross section that scales linearly with energy and where
// the differential cross section is proportional to the out-going
// neutrino energy.

class LinearCrossSections : public NeutrinoCrossSections {
  private:
    const squids::Const units;
    const double GF = 1.16639e-23; // eV^-2
    const double mp = 938.272e6; // proton mass eV
    const double CC_to_NC; // proportion of CC to NC which goes from 0 to 1.
  public :
    LinearCrossSections(double CC_to_NC):CC_to_NC(CC_to_NC){assert( CC_to_NC <= 1.0  && CC_to_NC >= 0.0 );}
    LinearCrossSections():LinearCrossSections(0.5){}
    double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
    double SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
  };
} // close nusquids namespace

#endif
