#include "exCross.h"

namespace nusquids{
  double LinearCrossSections::TotalCrossSection(double Enu, NeutrinoFlavor flavor,
								       NeutrinoType neutype, Current current) const {
    if(current == NeutrinoCrossSections::Current::CC)
      return 1.e-7*2.0*(1.-CC_to_NC)*GF*GF*Enu*mp;
    else if(current == NeutrinoCrossSections::Current::NC)
      return 1.e-7*2.0*CC_to_NC*GF*GF*Enu*mp;
    else
      throw std::runtime_error("Invalid current");
  }

  double LinearCrossSections::SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const {
    if (E2 > E1)
      return 0.0;
    else if (E2 < 0)
      return 0.0;
    else
      return (TotalCrossSection(E1,flavor,neutype,current)/(E1/units.GeV))*(2.*E2/E1);
  }
} // close nusquids namespace

