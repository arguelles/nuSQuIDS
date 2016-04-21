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
 
class NeutrinoDISCrossSectionsFromTablesExtended : public NeutrinoCrossSections {
    private :
      /// \brief True if the class has being initialized
      bool is_init = false;
      /// \brief Minimum neutrino energy.
      double Emin;
      /// \brief Maximum neutrino energy.
      double Emax;
      /// \brief Number of divisions.
      unsigned int div;
      /// \brief GeV in eV
      const double GeV = 1.0e9;

      /// \brief Stores the neutrino charge current differential cross section.
      marray<double,4> dsde_CC_data;
      /// \brief Stores the neutrino neutral current differential cross section.
      marray<double,4> dsde_NC_data;
      /// \brief Stores the array of the log energies of the data tables.
      std::vector<double> logE_data_range;

      /// \brief GSL interpolator for the total cross section.
      marray<gsl_spline *,3> xs_inter;
      /// \brief GSL interpolator accelerator for the total cross section.
      marray<gsl_interp_accel *,3> xs_acc;

      /// \brief Bilinear interpolator
      /// \details Used by DifferentialCrossSectionl() to interpolate the differential cross section.
      double LinInter(double,double,double,double,double) const;
      /// \brief Null Double Differential Cross section
      virtual double DoubleDifferentialCrossSection(double E, double x, double y, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override
      {
        return 0.;
      }
    public :
      /// \brief Default constructor
      // NeutrinoCrossSections(){};
      /// \brief Detauls destructor
      virtual ~NeutrinoDISCrossSectionsFromTablesExtended();
      /// \brief Constructor for a given energy range
      /// \details Calcualte all relevant cross section in the nodes setting a logarithmic scale
      NeutrinoDISCrossSectionsFromTablesExtended(){Init();}
      /// \brief Initializer for a given energy range
      /// \details Calcualte all relevant cross section in the nodes setting a logarithmic scale
      void Init();

      /// \brief Returns the total neutrino cross section
      /// \details Used to interpolate the total cross sections.
      double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
      /// \brief Returns the Differential cross section with respect to the outgoing lepton energy.
      /// \param E1 Incident lepton energy.
      /// \param E2 Outgoing lepton energy.
      /// \param flavor Flavor index.
      /// \param neutype Can be either neutrino or antineutrino.
      /// \param current Can be either CC or NC.
      /// \return The cross section in cm^2 GeV^-1.
      double SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
      /// \brief Returns the number of energy nodes.
      unsigned int GetNumE() const {return div;}
      /// \brief Returns the minimum energy in [eV]
      double GetEmin() const {return Emin;}
      /// \brief Returns the maximum energy in [eV]
      double GetEmax() const {return Emax;}
      /// \brief Returns true if the object is initialized.
      bool IsInit() const {return is_init;}
};


  

}



#endif
