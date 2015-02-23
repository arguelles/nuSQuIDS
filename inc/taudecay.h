#ifndef __TAUDECAY_H
#define __TAUDECAY_H

#include "tools.h"
#include <string>
#include <cmath>
#include <math.h>

/*
 * FORMULAES implemented out of
 *
 * Tau neutrinos underground
 * PHYSICAL REVIEW D 62 123001
 * SHARADA IYER DUTTA, MARY HALL RENO, AND INA SARCEVIC
 *
 */

// Tau decay spectra

namespace nusquids{

/// \class TauDecaySpectra
/// \brief Implements all formulaes necesary to tau decay kinematics.
/// \details For details on the formulaes read "Tau neutrinos underground,
/// PRD 62 123001, Hall Reno et al.
class TauDecaySpectra{
  private:
    /// \brief Stores the tau polarization which is set by SetParameters()
    /// @see SetParameters
    double TauPolarization;

    /// \brief Ratio of pion mass to tau mass squared.
    /// @see SetParameters
    double RPion;
    /// \brief Ratio of rho mass to tau mass squared.
    /// @see SetParameters
    double RRho;
    /// \brief Ratio of A1 meson mass to tau mass squared.
    /// @see SetParameters
    double RA1;

    /// \brief Tau decay to leptons branching ratio.
    /// @see SetParameters
    double BrLepton;
    /// \brief Tau decay to hadrons branching ratio.
    /// @see SetParameters
    double BrHadron;
    /// \brief Tau decay to pion branching ratio.
    /// @see SetParameters
    double BrPion;
    /// \brief Tau decay to rho branching ratio.
    /// @see SetParameters
    double BrRho;
    /// \brief Tau decay to A1 branching ratio.
    /// @see SetParameters
    double BrRA1;

    /// \brief Sets mass ratios, branching ratios, and polarizations.
    void SetParameters();
  protected:
    /// \brief Calculates the differential spectrum for tau to leptons with respect to z = E_nu/E_tau.
    /// @param E_tau Tau energy.
    /// @param E_nu Tau neutrino energy.
    double TauDecayToLepton(double E_tau,double E_nu) const;
    /// \brief Calculates the differential spectrum for tau to pions with respect to z = E_nu/E_tau.
    /// @param E_tau Tau energy.
    /// @param E_nu Tau neutrino energy.
    double TauDecayToPion(double E_tau,double E_nu) const;
    /// \brief Calculates the differential spectrum for tau to rhos with respect to z = E_nu/E_tau.
    /// @param E_tau Tau energy.
    /// @param E_nu Tau neutrino energy.
    double TauDecayToRho(double E_tau,double E_nu) const;
    /// \brief Calculates the differential spectrum for tau to A1 with respect to z = E_nu/E_tau.
    /// @param E_tau Tau energy.
    /// @param E_nu Tau neutrino energy.
    double TauDecayToA1(double E_tau,double E_nu) const;
    /// \brief Calculates the differential spectrum for tau to hadrons with respect to z = E_nu/E_tau.
    /// @param E_tau Tau energy.
    /// @param E_nu Tau neutrino energy.
    double TauDecayToHadron(double E_tau,double E_nu) const;
    /// \brief Calculates the differential spectrum for tau to all channels with respect to z = E_nu/E_tau.
    /// @param E_tau Tau energy.
    /// @param E_nu Tau neutrino energy.
    double TauDecayToAll(double E_tau,double E_nu) const;
  private:
    /// \brief Stores the differential spectrum with respect to the incoming neutrino
    /// energy for all channels.
    Table dNdEnu_All_tbl;
    /// \brief Stores the differential spectrum with respect to the incoming neutrino
    /// energy for leptonic channels.
    Table dNdEnu_Lep_tbl;

    /// \brief Stores the differential spectrum with respect to the outgoing lepton
    /// energy for all channels.
    Table dNdEle_All_tbl;
    /// \brief Stores the differential spectrum with respect to the outgoing lepton
    /// energy for leptonic channels.
    Table dNdEle_Lep_tbl;

    /// \brief Minimum energy in the array.
    double Emin;
    /// \brief Maximum energy in the array.
    double Emax;
    /// \brief Number of energy divisions.
    unsigned int div;
  public :
    /// \brief Detault empty constructor.
    TauDecaySpectra();
    /// \brief Constructor for a given energy range.
    /// @param emin Minimum energy [eV]
    /// @param emax Maximum energy [eV]
    /// @param div Number of divisions
    /// \details Construct the tables from \c emin to \c emax with
    /// number of divisions \c div using a logarithmic scale.
    TauDecaySpectra(double emin,double emax,unsigned int div);
    /// \brief Initializer for a given energy range.
    /// @param emin Minimum energy [eV]
    /// @param emax Maximum energy [eV]
    /// @param div Number of divisions
    /// \details Construct the tables from \c emin to \c emax with
    /// number of divisions \c div using a logarithmic scale.
    void Init(double emin,double emax,unsigned int div);

    /// \brief Returns the differential spectrum with respect to the incoming neutrino energy for
    /// all decay channels. Returned in units of GeV^-1
    /// @param i_enu Initial energy node index.
    /// @param i_ele Outgoing energy node index.
    double dNdEnu_All(unsigned int i_enu,unsigned int i_ele) const;
    /// \brief Returns the differential spectrum with respect to the incoming neutrino energy for
    /// leptonic decay channel. Returned in units of GeV^-1.
    /// @param i_enu Initial energy node index.
    /// @param i_ele Outgoing energy node index.
    double dNdEnu_Lep(unsigned int i_enu,unsigned int i_ele) const;

    /// \brief Returns the differential spectrum with respect to the outgoing lepton energy for
    /// all decay channels. Returned in units of GeV^-1
    /// @param i_enu Initial energy node index.
    /// @param i_ele Outgoing energy node index.
    double dNdEle_All(unsigned int i_enu,unsigned int i_ele) const;
    /// \brief Returns the differential spectrum with respect to the outgoing neutrino energy for
    /// leptonic decay channel. Returned in units of GeV^-1
    /// @param i_enu Initial energy node index.
    /// @param i_ele Outgoing energy node index.
    double dNdEle_Lep(unsigned int i_enu,unsigned int i_ele) const;
};

}// close namespace

#endif
