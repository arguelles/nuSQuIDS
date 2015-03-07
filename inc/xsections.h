#ifndef __XSECTIONS_H
#define __XSECTIONS_H

#include "tools.h"
#include "marray.h"
#include "global.h"
#include <string>
#include <cmath>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <stdexcept>

// neutrino cross sections

namespace nusquids{

/// \class NeutrinoCrossSections
/// \brief Tabulates and interpolates all cross sections for a given energy array.
/// \details The cross section tables are supplied data/xsections/ and bilinear
/// interpolation is performed on the logarithm of the energy.
class NeutrinoCrossSections{
    public:
      /// \brief Neutrino flavors.
      enum NeutrinoFlavor {electron = 0, muon = 1, tau = 2, sterile = 4};
      /// \brief Neutrino types.
      enum NeutrinoType {neutrino = 0, antineutrino = 1};
    private :
      /// \brief True if the class has being initialized
      bool is_init = false;
      /// \brief Minimum neutrino energy.
      double Emin;
      /// \brief Maximum neutrino energy.
      double Emax;
      /// \brief Number of divisions.
      unsigned int div;
      /// \brief Stores the neutrino charge current differential cross section evaluated at the nodes.
      marray<double,4> dsde_CC_tbl;
      /// \brief Stores the neutrino neutral current differential cross section evaluated at the nodes.
      marray<double,4> dsde_NC_tbl;
      /// \brief Stores the total neutrino charge current cross section evaluated at the nodes.
      marray<double,3> sigma_CC_tbl;
      /// \brief Stores the total neutrino neutral current cross section evaluated at the nodes.
      marray<double,3> sigma_NC_tbl;

      /// \brief Stores the neutrino charge current differential cross section.
      marray<double,4> dsde_CC_data;
      /// \brief Stores the neutrino neutral current differential cross section.
      marray<double,4> dsde_NC_data;
      /// \brief Stores the total neutrino charge current cross section.
      //marray<double,3> sigma_CC_data;
      /// \brief Stores the total neutrino neutral current cross section.
      //marray<double,3> sigma_NC_data;
      /// \brief Stores the array of the log energyes of the data tables.
      std::vector<double> logE_data_range;

      /// \brief GSL interpolator for the total cross section.
      marray<gsl_spline *,3> xs_inter;
      /// \brief GSL interpolator accelerator for the total cross section.
      marray<gsl_interp_accel *,3> xs_acc;

      /// \brief Bilinear interpolator
      /// \details Used by DifferentialCrossSectionl() to interpolate the differential cross section.
      double LinInter(double,double,double,double,double) const;
    protected:
      /// \brief Interaction current type
      enum Current { CC, NC };
      /// \brief Returns the total neutrino cross section
      /// \details Used to interpolate the total cross sections.
      virtual double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const;
      /// \brief Returns the Differential cross section with respect to the outgoing lepton energy.
      /// \details The cross section will be returned in cm^2 GeV^-1.
      /// @param E1 Incident lepton energy.
      /// @param E2 Outgoing lepton energy.
      /// @param flavor Flavor index.
      /// @param neutype Can be either neutrino or antineutrino.
      /// @param current Can be either CC or NC.
      virtual double DifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current)const;
    public :
      /// \brief Default constructor
      // NeutrinoCrossSections(){};
      /// \brief Detauls destructor
      virtual ~NeutrinoCrossSections();
      /// \brief Constructor for a given energy range
      /// @param emin Minimum neutrino energy in eV.
      /// @param emax Maximum neutrino energy in eV.
      /// @param div Number of divisions in the logarithmic scale.
      /// \details Calcualte all relevant cross section in the nodes setting a logarithmic scale
      NeutrinoCrossSections(double emin,double emax,unsigned int div);
      /// \brief Initializer for a given energy range
      /// @param emin Minimum neutrino energy in eV.
      /// @param emax Maximum neutrino energy in eV.
      /// @param div Number of divisions in the logarithmic scale.
      /// \details Calcualte all relevant cross section in the nodes setting a logarithmic scale
      void Init(double emin,double emax,unsigned int div);

      /// \brief Returns the diferencial charge current cross section.
      /// @param i_enu Index of the incoming neutrino energy node.
      /// @param i_ele Index of the outgoing lepton energy node.
      /// @param flv Flavor of the neutrino: 0:electron, 1:muon, 2: tau.
      /// @param neutype 0:neutrino, 1:antineutrino.
      double dsde_CC(unsigned int i_enu, unsigned int i_ele, NeutrinoFlavor flv, NeutrinoType neutype) const;
      /// \brief Returns the diferencial neutral cross section.
      /// @param i_enu Index of the incoming neutrino energy node.
      /// @param i_ele Index of the outgoing lepton energy node.
      /// @param flv Flavor of the neutrino: 0:electron, 1:muon, 2: tau.
      /// @param neutype 0:neutrino, 1:antineutrino.
      double dsde_NC(unsigned int i_enu, unsigned int i_ele, NeutrinoFlavor flv, NeutrinoType neutype) const;
      /// \brief Returns the total charge cross section.
      /// @param i_enu Index of the incoming neutrino energy node.
      /// @param flv Flavor of the neutrino: 0:electron, 1:muon, 2: tau.
      /// @param neutype 0:neutrino, 1:antineutrino.
      double sigma_CC(unsigned int i_enu,NeutrinoFlavor flv,NeutrinoType neutype) const;
      /// \brief Returns the total neutral cross section.
      /// @param i_enu Index of the incoming neutrino energy node.
      /// @param flv Flavor of the neutrino: 0:electron, 1:muon, 2: tau.
      /// @param neutype 0:neutrino, 1:antineutrino.
      double sigma_NC(unsigned int i_enu,NeutrinoFlavor flv,NeutrinoType neutype) const;

      /// \brief Returns the number of energy nodes.
      unsigned int GetNumE() const {return div + 1;}
      /// \brief Returns the minimum energy in [eV]
      double GetEmin() const {return Emin;}
      /// \brief Returns the maximum energy in [eV]
      double GetEmax() const {return Emax;}
      /// \brief Returns true if the object is initialized.
      bool IsInit() const {return is_init;}
};

} // close namespace

#endif
