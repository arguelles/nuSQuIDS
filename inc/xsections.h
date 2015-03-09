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
    enum NeutrinoFlavor {electron = 0, muon = 1, tau = 2, sterile = 4};
    /// \brief Neutrino types.
    enum NeutrinoType {neutrino = 0, antineutrino = 1};
    /// \brief Interaction current type
    enum Current { CC, NC };
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
    virtual double DifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const = 0;
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

      /*
      /// \brief Stores the neutrino charge current differential cross section evaluated at the nodes.
      marray<double,4> dsde_CC_tbl;
      /// \brief Stores the neutrino neutral current differential cross section evaluated at the nodes.
      marray<double,4> dsde_NC_tbl;
      /// \brief Stores the total neutrino charge current cross section evaluated at the nodes.
      marray<double,3> sigma_CC_tbl;
      /// \brief Stores the total neutrino neutral current cross section evaluated at the nodes.
      marray<double,3> sigma_NC_tbl;
      */

      /// \brief Stores the neutrino charge current differential cross section.
      marray<double,4> dsde_CC_data;
      /// \brief Stores the neutrino neutral current differential cross section.
      marray<double,4> dsde_NC_data;
      /// \brief Stores the array of the log energyes of the data tables.
      std::vector<double> logE_data_range;

      /// \brief GSL interpolator for the total cross section.
      marray<gsl_spline *,3> xs_inter;
      /// \brief GSL interpolator accelerator for the total cross section.
      marray<gsl_interp_accel *,3> xs_acc;

      /// \brief Bilinear interpolator
      /// \details Used by DifferentialCrossSectionl() to interpolate the differential cross section.
      double LinInter(double,double,double,double,double) const;
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
      double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const;
      /// \brief Returns the Differential cross section with respect to the outgoing lepton energy.
      /// \details The cross section will be returned in cm^2 GeV^-1.
      /// @param E1 Incident lepton energy.
      /// @param E2 Outgoing lepton energy.
      /// @param flavor Flavor index.
      /// @param neutype Can be either neutrino or antineutrino.
      /// @param current Can be either CC or NC.
      double DifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const;

      /*
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
      */

      /// \brief Returns the number of energy nodes.
      unsigned int GetNumE() const {return div;}
      /// \brief Returns the minimum energy in [eV]
      double GetEmin() const {return Emin;}
      /// \brief Returns the maximum energy in [eV]
      double GetEmax() const {return Emax;}
      /// \brief Returns true if the object is initialized.
      bool IsInit() const {return is_init;}
};


} // close namespace

#endif
