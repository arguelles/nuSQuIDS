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
//#include <string>

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>


namespace nusquids{


// NeutrinoDISCrossSectionsFromTablesExtended
// Is basically a copy of the default NeutrinoDISCrossSectionsFromTables but instead of 
// returning error when the energy is lower that the low energy value in the tables it returns zero.
// This is effectively true for the range of energies given in the nuSQuIDS default tables(Emin=1e2GeV)
// This allows to use a wider range of the energy to compute atmospheric oscillations as it's shown in the 
// main.cpp example.
  
  class NeutrinoDISCrossSectionsFromTablesExtended : public NeutrinoDISCrossSectionsFromTables {
  public :
    NeutrinoDISCrossSectionsFromTablesExtended():NeutrinoDISCrossSectionsFromTables(){}
    double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
    double SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
  };
  
  
  
  
}



#endif
