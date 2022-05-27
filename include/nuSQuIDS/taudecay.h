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


#ifndef NUSQUIDS_TAUDECAY_H
#define NUSQUIDS_TAUDECAY_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>

#include "nuSQuIDS/marray.h"
#include "nuSQuIDS/aligned_alloc.h"

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
    ///aligned arrays should be aligned to at least this many items
    constexpr static size_t preferred_alignment=4;
    constexpr static uint8_t log2(size_t n) {
      return((n>1) ? log2(n/2)+1 : 0);
    }
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
  public:
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
    marray<double,2,aligned_allocator<double>> dNdEnu_All_tbl;
    /// \brief Stores the differential spectrum with respect to the incoming neutrino
    /// energy for leptonic channels.
    marray<double,2,aligned_allocator<double>> dNdEnu_Lep_tbl;

    /// \brief Stores the differential spectrum with respect to the outgoing lepton
    /// energy for all channels.
    marray<double,2,aligned_allocator<double>> dNdEle_All_tbl;
    /// \brief Stores the differential spectrum with respect to the outgoing lepton
    /// energy for leptonic channels.
    marray<double,2,aligned_allocator<double>> dNdEle_Lep_tbl;

  public :
    /// \brief Detault empty constructor.
    TauDecaySpectra():
      dNdEnu_All_tbl(aligned_allocator<double>{log2(preferred_alignment*sizeof(double))}),
      dNdEnu_Lep_tbl(aligned_allocator<double>{log2(preferred_alignment*sizeof(double))}),
      dNdEle_All_tbl(aligned_allocator<double>{log2(preferred_alignment*sizeof(double))}),
      dNdEle_Lep_tbl(aligned_allocator<double>{log2(preferred_alignment*sizeof(double))})
      {
        SetParameters();
      }
    /// \brief Constructor for a given energy range.
    /// @param E_range Energy nodes where the cross section will be calculated. [eV]
    /// \details Construct the tables on a rectangular grid given by E_range X E_range.
    TauDecaySpectra(marray<double,1> E_range);
    /// \brief Initializer for a given energy range.
    /// @param E_range Energy nodes where the cross section will be calculated. [eV]
    /// \details Construct the tables on a rectangular grid given by E_range X E_range.
    void Init(marray<double,1> E_range);

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
    /// \brief Returns tau to lepton branching ratio.
    double GetTauToLeptonBranchingRatio() const{
      return BrLepton;
    }
    double GetTauToHadronBranchingRatio() const{
      return BrHadron;
    }

    double DecayLenght(double Etau) const { return 87e-6*Etau/1.77686; }; //decay lenght of tau is 87 micrometers



};

}// close namespace

#endif
