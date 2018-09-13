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


#ifndef NUSQUIDS_XSECTIONS_H
#define NUSQUIDS_XSECTIONS_H

#include "version.h"
#include "tools.h"
#include "marray.h"
#include "global.h"
#include <string>
#include <cmath>
#include <stdexcept>
#include "SQuIDS/const.h"

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
    /// \brief Returns the double differential cross section with respect to x and y.
    /// \details The cross section will be returned in cm^2 GeV^-1. As this cross sections is not requiered and thus not called by the program
    /// its implementation is optional.
    /// @param E Incident lepton energy.
    /// @param x bjorken-x.
    /// @param y bjorken-y.
    /// @param flavor Flavor index.
    /// @param neutype Can be either neutrino or antineutrino.
    /// @param current Can be either CC or NC.
    virtual double DoubleDifferentialCrossSection(double E, double x, double y, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const {
      throw std::runtime_error("NeutrinoCrossSections::Error::DoubleDifferentialCrossSection is not implemented.");
      return 0;
    }
};

/// \class NullCrossSections
/// \brief Simple class that defines a null cross section.
class NullCrossSections: public NeutrinoCrossSections {
  public :
    NullCrossSections(){}
    /// \brief Returns the total neutrino cross section
    double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override { return 0;}
    /// \brief Returns the Differential cross section with respect to the outgoing lepton energy.
    double SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override { return 0;}
    /// \brief Returns the double differential cross section with respect to x and y.
    double DoubleDifferentialCrossSection(double E, double x, double y, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override { return 0;}
};

/// \class NeutrinoDISCrossSectionsFromTables
/// \brief Tabulates and interpolates all cross sections for a given energy array.
/// \details The cross section tables are supplied data/xsections/ and bilinear
/// interpolation is performed on the logarithm of the energy.
class NeutrinoDISCrossSectionsFromTables : public NeutrinoCrossSections {
    protected :
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

      /// \brief The neutrino charged current total cross section
      ///
      ///Indices are neutrino/anti-neutrino, flavor, energy
      marray<double,3> s_CC_data;
      /// \brief The neutrino neutral current total cross section
      ///
      ///Indices are neutrino/anti-neutrino, flavor, energy
      marray<double,3> s_NC_data;
      /// \brief The neutrino charged current differential cross section.
      ///
      ///Indices are neutrino/anti-neutrino, flavor, incident energy, out-going energy
      marray<double,4> dsde_CC_data;
      /// \brief The neutrino neutral current differential cross section.
      ///
      ///Indices are neutrino/anti-neutrino, flavor, incident energy, out-going energy
      marray<double,4> dsde_NC_data;
      /// \brief Stores the array of the log energies of the data tables.
      std::vector<double> logE_data_range;

      /// \brief Linear interpolator
      /// \details Used by DifferentialCrossSectionl() to interpolate the differential cross section.
      /// \param x abscissa value at which to compute the interpolated function value
      /// \param xM closest tabulated abscissa value smaller than x
      /// \param xP closest tabulated abscissa value larger than x
      /// \param yM tablated ordinate value at xM
      /// \param yP tablated ordinate value at xP
      double LinInter(double x,double xM,double xP,double yM,double yP) const;
      /// \brief Null Double Differential Cross section
      virtual double DoubleDifferentialCrossSection(double E, double x, double y, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override
      {
        return 0.;
      }
      /// \brief Reads cross section data from a collection of text files
      /// \param root the base file path
      void ReadText(std::string root);
    public :
      virtual ~NeutrinoDISCrossSectionsFromTables();
      /// \brief Default construct with built-in tables
      NeutrinoDISCrossSectionsFromTables();
      /// \brief Construct with data from the given file(s)
      /// \details If the specified path is a file it will be assumed to be an
      ///          HDF5 file, otherwise it must be a prefix of several text files
      ///          with related names which contain the cross sections. All cross
      ///          sections should be in units of cm^2 (or cm^2/GeV). 
      ///
      ///          If the input is an HDF5 file it must contain the following:
      ///          - Two attributes `Emin` and `Emax` which specify the minimum
      ///            and maximum tabulated energies, in eV, respectively. 
      ///          - Two 3 dimensional datasets named `s_CC` and `s_NC` which 
      ///            contain the total cross sections for all particle types, 
      ///            flavors, and energies for charged current and neutral 
      ///            current interactions, respectively. The first dimension of
      ///            each dataset is indexed by particle type (neutrino or 
      ///            antineutrino) and must have extent 2. The second dimension
      ///            is indexed by flavor (0 = electron, 1 = muon, 2 = tau) and
      ///            must have extent 3. The final dimension is indexed by 
      ///            incoming neutrino energy and must have an extent of at 
      ///            least 2. The energiesfor which the cross section is 
      ///            tabulated must be logarthmically spaced. 
      ///          - Two 4 dimensional tables named `dsDE_CC` and `dsDE_NC`
      ///            which contain the differential cross sections in out-going
      ///            lepton energy for charged current and neutral interactions, 
      ///            respectively. The first three dimensions are the same as 
      ///            for the total cross section tables. The final dimension is
      ///            indexed by the out-going lepton energy and must have the 
      ///            same extent and correspond to the same energy values as the
      ///            incoming neutrino dimension. 
      ///
      ///          If the input is not a single HDF5 file is must be a set of 
      ///          text files whose names are the same except that they have the 
      ///          suffixes:
      ///          - 'sigma_CC.dat' for the total charged current cross section 
      ///          - 'sigma_NC.dat' for the total neutral current cross section 
      ///          - 'dsde_CC.dat' for the differential charged current cross section 
      ///          - 'dsde_NC.dat' for the differential neutral current cross section 
      ///
      ///          The total cross section files must contain on each line first
      ///          the energy _in GeV_ at which the cross sections are evaluated,
      ///          then pairs of cross section values for all three flavors in 
      ///          the order electron, muon, tau, with each pair consisting of
      ///          the cross sections for neutrinos and antineutrinos. The energies
      ///          must be the same in all files, must be listed in ascending 
      ///          order, and must be logarithmically spaced. 
      ///
      ///          The differential cross section files must contain the data
      ///          as the total cross section files, except that a column is 
      ///          inserted after the first for the out-going lepton energy. The
      ///          set of out-going energies must be the same as the incoming 
      ///          energies, and the lines in the files must be lexicographically
      ///          ordered by the two energies. 
      /// \param path Either a path to a single HDF5 file which contains all
      ///             cross section data, or the base path for a set of text
      ///             files containing the various cross sections
      NeutrinoDISCrossSectionsFromTables(std::string path);

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
      /// \brief Write the cross sections to an HDF5 file
      /// \param path the path to which the file should be written
      void WriteHDF(std::string path) const;
      /// \brief Write the cross sections to a set of text files
      /// \param basePath the base path name for the output files;
      ///                 suffixes will be automatically appended.
      void WriteText(std::string basePath) const;
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
  const squids::Const constants;
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
