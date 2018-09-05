#include "exCross.h"

namespace nusquids{
  double LinearCrossSections::TotalCrossSection(double Enu, NeutrinoFlavor flavor,
								       NeutrinoType neutype, Current current) const {
    return 1.e-7*GF*GF*Enu*mp;
  }

  double LinearCrossSections::SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const {
    if (E2 > E1)
      return 0.0;
    else if (E2 < 0)
      return 0.0;
    else
      return (TotalCrossSection(E1,flavor,neutype,current)/E1)*(2.*E2/E1);
  }
} // close nusquids namespace

