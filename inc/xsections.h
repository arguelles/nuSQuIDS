#ifndef __XSECTIONS_H
#define __XSECTIONS_H

#include "tools.h"
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
    private :
      /// \brief Interaction current type
      enum Current { CC, NC };
      /// \brief Minimum neutrino energy.
      double Emin;
      /// \brief Maximum neutrino energy.
      double Emax;
      /// \brief Number of divisions.
      int div;

      /// \brief Stores the neutrino charge current differential cross section evaluated at the nodes.
      Table dsde_CC_tbl;
      /// \brief Stores the neutrino neutral current differential cross section evaluated at the nodes.
      Table dsde_NC_tbl;
      /// \brief Stores the total neutrino charge current cross section evaluated at the nodes.
      Table sigma_CC_tbl;
      /// \brief Stores the total neutrino neutral current cross section evaluated at the nodes.
      Table sigma_NC_tbl;

      /// \brief Stores the neutrino charge current differential cross section.
      Table dsde_CC_data;
      /// \brief Stores the neutrino neutral current differential cross section.
      Table dsde_NC_data;
      /// \brief Stores the total neutrino charge current cross section.
      Table sigma_CC_data;
      /// \brief Stores the total neutrino neutral current cross section.
      Table sigma_NC_data;

      /// \brief Bilinear interpolator
      /// \details Used by DSDEInter() to interpolate the differential cross section.
      double LinInter(double,double,double,double,double);
      /// \brief Spline interpolator
      /// \details Used to interpolate the total cross sections.
      double SigmaInter(double, gsl_spline *, gsl_interp_accel *);
      /// \brief  interpolator
      /// \details
      /// @param E1 Incident lepton energy.
      /// @param E2 Outgoing lepton energy.
      /// @param flavor Flavor index.
      /// @param logE_range Array containing the logarithm of the energy nodes.
      /// @param current Can be either "CC" or "NC".
      double DSDEInter(double E1,double E2,int flavor ,std::vector<double> logE_range,Current current);
    public :
      /// \brief Default constructor
      NeutrinoCrossSections();
      /// \brief Constructor for a given energy range
      /// @param emin Minimum neutrino energy in eV.
      /// @param emax Maximum neutrino energy in eV.
      /// @param div Number of divisions in the logarithmic scale.
      /// \details Calcualte all relevant cross section in the nodes setting a logarithmic scale
      NeutrinoCrossSections(double emin,double emax,int div);
      /// \brief Initializer for a given energy range
      /// @param emin Minimum neutrino energy in eV.
      /// @param emax Maximum neutrino energy in eV.
      /// @param div Number of divisions in the logarithmic scale.
      /// \details Calcualte all relevant cross section in the nodes setting a logarithmic scale
      void Init(double emin,double emax,int div);

      /// \brief Returns the diferencial charge current cross section.
      /// @param i_enu Index of the incoming neutrino energy node.
      /// @param i_ele Index of the outgoing lepton energy node.
      /// @param flv Flavor of the neutrino: 0:electron, 1:muon, 2: tau.
      /// @param neutype 0:neutrino, 1:antineutrino.
      double dsde_CC(int i_enu, int i_ele, int flv, int neutype);
      /// \brief Returns the diferencial neutral cross section.
      /// @param i_enu Index of the incoming neutrino energy node.
      /// @param i_ele Index of the outgoing lepton energy node.
      /// @param flv Flavor of the neutrino: 0:electron, 1:muon, 2: tau.
      /// @param neutype 0:neutrino, 1:antineutrino.
      double dsde_NC(int i_enu, int i_ele, int flv, int neutype);
      /// \brief Returns the total charge cross section.
      /// @param i_enu Index of the incoming neutrino energy node.
      /// @param flv Flavor of the neutrino: 0:electron, 1:muon, 2: tau.
      /// @param neutype 0:neutrino, 1:antineutrino.
      double sigma_CC(int i_enu,int flv,int neutype);
      /// \brief Returns the total neutral cross section.
      /// @param i_enu Index of the incoming neutrino energy node.
      /// @param flv Flavor of the neutrino: 0:electron, 1:muon, 2: tau.
      /// @param neutype 0:neutrino, 1:antineutrino.
      double sigma_NC(int i_enu,int flv,int neutype);
};

} // close namespace

#endif
