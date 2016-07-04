#include "exCross.h"


namespace nusquids{


  double NeutrinoDISCrossSectionsFromTablesExtended::TotalCrossSection(double Enu, NeutrinoFlavor flavor,
								       NeutrinoType neutype, Current current) const{
    // we assume that sterile neutrinos are trully sterile
    if (not (flavor == electron or flavor == muon or flavor == tau))
      return 0.0;
    if (Enu < Emin)
      return std::numeric_limits<double>::min();
    else
      return NeutrinoDISCrossSectionsFromTables::TotalCrossSection(Enu, flavor, neutype, current);
  }
  
  double NeutrinoDISCrossSectionsFromTablesExtended::SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const{
    // we assume that sterile neutrinos are trully sterile
    if (not (flavor == electron or flavor == muon or flavor == tau))
      return 0.0;
    if (E1 < Emin || E2<Emin)
      return std::numeric_limits<double>::min();
    else
      return NeutrinoDISCrossSectionsFromTables::SingleDifferentialCrossSection(E1, E2, flavor, neutype, current);
  }
    
}

