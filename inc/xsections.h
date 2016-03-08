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


#ifndef __XSECTIONS_H
#define __XSECTIONS_H

#include "version.h"
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
/// \brief General neutrino cross section class.
class NeutrinoCrossSections{
  public:
    /// \brief Neutrino flavors.
    enum NeutrinoFlavor {electron = 0, muon = 1, tau = 2, sterile = 3};
    /// \brief Neutrino types.
    enum NeutrinoType {neutrino = 0, antineutrino = 1};
    /// \brief Interaction current type
    enum Current { CC, NC, GR };
    /// \brief Returns the total neutrino cross section
    /// \details Used to interpolate the total cross sections.
    virtual double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const = 0;
    /// \brief Returns the Differential cross section with respect to the outgoing lepton energy.
    /// \details The cross section will be returned in cm^2 GeV^-1.
    /// @param E1 Incident lepton energy.
    /// @param E2 Outgoing lepton energy.
    /// @param flavor Flavor index.
    /// @param neutype Can be either neutrino or antineutrino.
    /// @param current Can be either CC or NC.
    virtual double SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const = 0;
    /// \brief Returns the Differential cross section with respect to the outgoing lepton energy.
    /// \details The cross section will be returned in cm^2 GeV^-1.
    /// @param E Incident lepton energy.
    /// @param x bjorken-x.
    /// @param y bjorken-y.
    /// @param flavor Flavor index.
    /// @param neutype Can be either neutrino or antineutrino.
    /// @param current Can be either CC or NC.
    virtual double DoubleDifferentialCrossSection(double E, double x, double y, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const = 0;
};

/// \class NeutrinoDISCrossSectionsFromTables
/// \brief Tabulates and interpolates all cross sections for a given energy array.
/// \details The cross section tables are supplied data/xsections/ and bilinear
/// interpolation is performed on the logarithm of the energy.
class NeutrinoDISCrossSectionsFromTables : public NeutrinoCrossSections {
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
      virtual ~NeutrinoDISCrossSectionsFromTables();
      /// \brief Constructor for a given energy range
      /// \details Calcualte all relevant cross section in the nodes setting a logarithmic scale
      NeutrinoDISCrossSectionsFromTables(){Init();}
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

/// \class NeutrinoGRCrossSection
/// \brief Implements electron-antineutrino/electron scattering
class GlashowResonanceCrossSection : public NeutrinoCrossSections {
public:
  virtual ~GlashowResonanceCrossSection();

  GlashowResonanceCrossSection();

  /// \brief Returns the total neutrino cross section
  double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
  /// \brief Returns the Differential cross section with respect to the outgoing lepton energy.
  /// \param E1 Incident lepton energy.
  /// \param E2 Outgoing lepton energy.
  /// \param flavor Flavor index. Must be 0
  /// \param neutype Must be antineutrino.
  /// \param current Must be GR.
  /// \return The cross section in cm^2 GeV^-1.
  double SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
  
  /// \warning Not implemented; should not be used.
  double DoubleDifferentialCrossSection(double E, double x, double y, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
  
  /// \brief Returns the fraction of final states with a lepton of the specified flavor
  /// \details For muon flavor, this is the ratio of the integrated single-differential cross-section to the total cross-section
  static double WDecayBranchingFraction(NeutrinoFlavor flavor){
    switch(flavor){
      case electron: return(B_Electron);
      case muon: return(B_Muon);
      case tau: return(B_Tau);
      default: return(0);
    }
  }

private:
  double fermi_scale;
  ///W mass
  double M_W;
  ///W full decay width
  double W_total;
  
  ///W^+ -> e^+ + \nu_e branching ratio
  static double B_Electron;
  ///W^+ -> \mu^+ + \nu_\mu branching ratio
  static double B_Muon;
  ///W^+ -> \tau^+ + \nu_\tau branching ratio
  static double B_Tau;
  ///W^+ -> hadrons branching ratio
  static double B_Hadronic;
};

} // close namespace

#endif
