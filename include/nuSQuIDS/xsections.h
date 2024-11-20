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

#include <string>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <SQuIDS/const.h>
#include "nuSQuIDS/version.h"
#include "nuSQuIDS/tools.h"
#include "nuSQuIDS/marray.h"

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
    /// \param Enu Incident lepton energy
    /// \param flavor Incident lepton flavor
    /// \param neutype Incident lepton matter/anti-matter type
    /// \param current Interaction type
    /// \return The cross section in cm^2
    virtual double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const = 0;
    /// \brief Returns the differential cross section with respect to the outgoing lepton energy.
    /// \details The cross section will be returned in cm^2 GeV^-1.
    /// \param E1 Incident lepton energy.
    /// \param E2 Outgoing lepton energy.
    /// \param flavor Flavor index.
    /// \param neutype Can be either neutrino or antineutrino.
    /// \param current Can be either CC or NC.
    virtual double SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const = 0;
    /// \brief Returns the double differential cross section with respect to x and y.
    /// \details The cross section will be returned in cm^2 GeV^-1. As this cross sections is not  
    /// generally used by the library, implementation of this function is optional. 
    /// \param E Incident lepton energy.
    /// \param x bjorken-x.
    /// \param y bjorken-y.
    /// \param flavor Flavor index.
    /// \param neutype Can be either neutrino or antineutrino.
    /// \param current Can be either CC or NC.
    virtual double DoubleDifferentialCrossSection(double E, double x, double y, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const {
      throw std::runtime_error("NeutrinoCrossSections::Error::DoubleDifferentialCrossSection is not implemented.");
      return 0;
    }
    
    /// \brief Returns the total neutrino cross section, averaged over the specified energy range.
    /// \param EMin Minimum incident lepton energy
    /// \param EMax Maximum incident lepton energy
    /// \param flavor Incident lepton flavor
    /// \param neutype Incident lepton matter/anti-matter type
    /// \param current Interaction type
    /// \return The average cross section from EMin to Emax in cm^2
    virtual double AverageTotalCrossSection(double EMin, double EMax, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const;
    /// \brief Returns the the differential cross section, with respect to the outgoing lepton energy, averaged over the specified energy range.
    /// \param E1 Incident lepton energy
    /// \param E2Min Minimum out-going lepton energy
    /// \param E2Max Maximum out-going lepton energy
    /// \param flavor Incident lepton flavor
    /// \param neutype Incident lepton matter/anti-matter type
    /// \param current Interaction type
    /// \return The average cross section from EMin to Emax in cm^2 GeV^-1
    virtual double AverageSingleDifferentialCrossSection(double E1, double E2Min, double E2Max, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const;
};

struct FlavorHash{
public:
  using argument_type=NeutrinoCrossSections::NeutrinoFlavor;
private:
  using underlying_type=typename std::underlying_type<argument_type>::type;
  using hash_type=std::hash<underlying_type>;
public:
  using result_type=typename hash_type::result_type;
  result_type operator()(argument_type arg) const{
    return hash_(arg);
  }
private:
  hash_type hash_;
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

/// \class NeutrinoDISCrossSectionsFromTables_V1
/// \brief Tabulates and interpolates all cross sections for a given energy array.
/// \details The cross section tables are supplied data/xsections/ and bilinear
/// interpolation is performed on the logarithm of the energy.
class NeutrinoDISCrossSectionsFromTables_V1 : public NeutrinoCrossSections {
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
      virtual ~NeutrinoDISCrossSectionsFromTables_V1();
      /// \brief Default construct with built-in tables
      NeutrinoDISCrossSectionsFromTables_V1();
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
      NeutrinoDISCrossSectionsFromTables_V1(std::string path);

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

/// \brief Interpolates tabulated DIS cross sections with cubic polynomials
/// \details Total cross sections are stored and interpolated in double-log
/// space, that is the base ten logarithm of the total cross section is 
/// tabulated in the base ten logarithm of the incoming neutrino energy. The
/// singly-differential cross section is uses the same axes, but its additional 
/// axis for out-going lepton energy is tabulated linearly in energy above the 
/// minimum tabulated incoming neutrino energy scaled by the difference in 
/// incoming neutrino energy and the minimum tabulated energy:
/// z = (E_out - E_min)/(E_in - E_min)
/// This alternate variable has a domain of [0,1], which has the useful property
/// that the physical boundary E_out <= E_in maps to z <= 1, (and E_out >= 0 
/// becomes z >= 0), so the full space of a rectangular table can be utilized, 
/// and interpolation does not need to address special cases of grid cells which 
/// have only some bounding points in the physical region. A linear mapping is
/// suitable because DIS cross sections have the helpful property of becoming 
/// smooth and even near linear for E_out << E_in, so there is little need for 
/// detailed tabulation around z = 0. Conversely, the largest contributions to 
/// the cross section occur near E_out = E_in - epsilon, which region is better 
/// sampled by a linear mapping than a logarithmic one. It is still advisable, 
/// when computing the singly-differential cross sections, if the maximum of the
/// cross section for a given E_in is between the last two samples in z to 
/// replace the sample at z = 1 (whose value should formally be zero) with the 
/// peak value. nuSQuiDS assumes implicitly that the cross section for 
/// E_out = E_in is zero, and so will not use it explicitly, and at high E_in 
/// where the cross section becomes sharply peaked, missing the peak between the 
/// last to samples in z can lead to a severe and systematic under-estimate of
/// the average cross section in that region, corresponding to an under-estimate
/// of absorption. 
class NeutrinoDISCrossSectionsFromTables : public NeutrinoCrossSections {
protected:
	/// \brief Minimum neutrino energy.
	double Emin;
	/// \brief Maximum neutrino energy.
	double Emax;
	///\brief Conversion factor from GeV to eV
	const double GeV = 1.0e9;
	
    ///\brief Cross section information for one neutrino flavor
    ///
    ///Assumed to use a common energy range.
    struct perFlavorData{
      ///Total cross section for charged-current neutrino interactions
      AkimaSpline s_CC_nu;
      ///Total cross section for neutral-current neutrino interactions
      AkimaSpline s_NC_nu;
      ///Total cross section for charged-current anti-neutrino interactions
      AkimaSpline s_CC_nubar;
      ///Total cross section for neutral-current anti-neutrino interactions
      AkimaSpline s_NC_nubar;
      
      ///Singly-differential cross section for charged-current neutrino interactions
      BiCubicInterpolator dsdy_CC_nu;
      ///Singly-differential cross section for neutral-current neutrino interactions
      BiCubicInterpolator dsdy_NC_nu;
      ///Singly-differential cross section for charged-current anti-neutrino interactions
      BiCubicInterpolator dsdy_CC_nubar;
      ///Singly-differential cross section for neutral-current anti-neutrino interactions
      BiCubicInterpolator dsdy_NC_nubar;
    };
  
    std::unordered_map<NeutrinoFlavor,std::shared_ptr<perFlavorData>,FlavorHash> xsData;
	
	///Read a total cross section table from whitespace-separated text
	///\param path the filesystem path from which to read input
	///\return a tuple of a spline interpolating the tabulated data, the minimum
	///        tabulated energy, and maximum tabulated energies. Energies are in
	///        natural units. 
	std::tuple<AkimaSpline,double,double> read1DInterpolationFromText(const std::string& path);
	
	///Read a singly-differential cross section table from whitespace-separated text
	///\param path the filesystem path from which to read input
	///\return a tuple of a bicubic interpolator for the tabulated data, the 
	///        minimum tabulated energy, and maximum tabulated energy.
	///        Energies are in natural units. 
	std::tuple<BiCubicInterpolator,double,double> read2DInterpolationFromText(const std::string& path);
  
    ///Read data for a single flavor from a set of text files with a common prefix
    ///\param prefix the common initial part of the text files' path and name
    ///\param erangeSrc a description of the source from which the energy range 
    ///                 has been determined, if any
    std::shared_ptr<perFlavorData> readFlavorText(const std::string& prefix, 
	                                              std::string& erangeSrc);
	
    ///Read data for a single flavor from an HDF5 group
    ///\param sourceLoc the group (or file) from which to read HDF datasets
    ///\param energies the common set of energy values used for all tables
    ///\param zs the common set of z values used for all tables
    std::shared_ptr<perFlavorData> readFlavorHDF5(hid_t sourceLoc, 
                                                  const marray<double,1>& energies, 
                                                  const marray<double,1>& zs);
	
    ///Write data for a single flavor to an HDF5 group
    ///\param destLoc the group (or file) to which to write HDF datasets
    ///\param data the data to write
    ///\param compressionLevel level of zlib compression to apply to datasets
    void writeFlavorHDF5(hid_t destLoc, const perFlavorData& data, unsigned int compressionLevel) const;
	
	///\brief Check whether a path refers to an HDF5 file
	bool isHDF(const std::string& path);
	
    ///\brief Get the data corresponding to a given flavor
    const perFlavorData& getFlavor(NeutrinoFlavor flavor) const;
	
public:
	///Construct a set of cross sections from a default set of tables
	NeutrinoDISCrossSectionsFromTables();

	///Construct a set of cross sections from either ASCII text tables or an 
	///HDF5 file. See readText and readHDF for descriptions of input formats. 
	///\param pathOrPrefix either common prefix for the paths of a family of 
	///                    text files or the path to a single HDF5 file
	NeutrinoDISCrossSectionsFromTables(std::string pathOrPrefix);
	
	double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
	
	double SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
	
	//Implementing this doesn't give any substantial benefit over the default implementation
	//double AverageTotalCrossSection(double EnuMin, double EnuMax, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
	
	double AverageSingleDifferentialCrossSection(double E1, double E2Min, double E2Max, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
	
    /// Read a set of cross sections from text files with a common path prefix. 
    ///
    /// This will fail messily and leave the object in an inconsistent state if 
    /// any file does not exist or has the wrong structure. 
    ///
    /// The input files must have consistent names baased on the common prefix:
    /// - ${prefix}${fl_id}nu_sigma_CC.dat : the total charged-current neutrino xs
    /// - ${prefix}${fl_id}nu_sigma_NC.dat : the total neutral-current neutrino xs
    /// - ${prefix}${fl_id}nubar_sigma_CC.dat : the total charged-current anti-neutrino xs
    /// - ${prefix}${fl_id}nubar_sigma_NC.dat : the total neutral-current anti-neutrino xs
    /// - ${prefix}${fl_id}nu_dsde_CC.dat : the singly-differential charged-current neutrino xs
    /// - ${prefix}${fl_id}nu_dsde_NC.dat : the singly-differential neutral-current neutrino xs
    /// - ${prefix}${fl_id}nubar_dsde_CC.dat : the singly-differential charged-current anti-neutrino xs
    /// - ${prefix}${fl_id}nubar_dsde_NC.dat : the singly-differential neutral-current anti-neutrino xs
    ///
    /// The fl_id component may be the empty string, in which case the data
    /// in the files is assumed to be valid for all active neutrino flavors. 
    /// Alternatively, there may be three sets of files with fl_id taking the
    /// values "electron_", "muon_", and "tau_", respectively, in which case the
    /// sets of files are interpreted as corresponding distinctly to the active 
    /// flavors in the obvious way. 
    ///
    /// Each file must cover the same range of incoming neutrino energies, with 
    /// at least three energies being sampled. 
    /// Each total cross section file must contain two whitespace-separated
    /// columns, incoming neutrino energy and cross section. Entries must be 
    /// sorted in ascending order of incoming neutrino energy. Energies must be
    /// expressed in GeV, and cross sections in square centimeters. 
    /// Each singly-differential cross section file must have three columns: 
    /// incoming neutrino energy, 'z', and cross section. Energies must be in
    /// GeV and cross sections in square centimeters, as before. 'z' is a 
    /// unitless variable which describes the out-going lepton energy. 
    /// 'z' is defined as:
    ///     z = (E_out-E_min)/(E_in-E_min)
    /// where E_min is the minimum tabulated incoming neutrino energy, E_in is
    /// the incoming neutrino energy for the table entry, and E_out is the 
    /// out-going lepton energy. Entries must be sorted in ascending order of
    /// incoming neutrino energy, and ascending order of 'z' for matching
    /// incoming energies. 'z' values must be uniformly sampled over the domain
    /// [0,1], and the same numbers of entries must be present for all incoming
    /// neutrino energies, i.e. the table must be rectangular. At least three 
    /// values of 'z' must be sampled. 'z' is not well-defined for the lowest 
    /// incoming neutrino energy in the table, E_in = E_min. This can be 
    /// ignored, as there can be no tabulated entry with E_out < E_in for this 
    /// energy (and thus a non-zero cross section), instead the same set of z 
    /// values as for other energies should be used, each with a cross section 
    /// of zero. 
    ///
    ///\param prefix the common filesystem path prefix for all input files
    void ReadText(const std::string& prefix);
    
    /// Write out a set of text files in the format expected by ReadText. 
    ///\param prefix the path prefix from which each output file's name should 
    ///              be derived
    void WriteText(const std::string& prefix) const;
	
    /// Read a set of cross sections from an HDF5 file. 
    ///
    /// This will fail messily and leave the object in an inconsistent state if 
    /// the file has the wrong structure. 
    ///
    /// The file must contain a set of ten datasets:
    /// - energies
    /// - zs
    /// - s_CC_nu
    /// - s_NC_nu
    /// - s_CC_nubar
    /// - s_NC_nubar
    /// - dsdy_CC_nu
    /// - dsdy_NC_nu
    /// - dsdy_CC_nubar
    /// - dsdy_NC_nubar
    ///
    /// Or, alternatively, the file may contain the energies and zs datasets
    /// directly and also contain three groups:
    /// - electron
    /// - muon
    /// - tau
    /// each containing its own set of the eight s_* and dsdy_* datasets.
    ///
    /// In the former case all total and singly differential cross sections will
    /// be assumed to be common to all active neutrino flavors, while in the
    /// latter case the cross sections will be distinctly assigned to the active 
    /// flavors in the obvious way. 
    ///
    /// energies and zs must be one dimensional datasets, containing the 
    /// base-ten logarithm of the energy values (in electronvolts) and the 
    /// values of the 'z' transformed variable used for the other eight tables.
    /// The z variable is defined as:
    ///     z = (E_out-E_min)/(E_in-E_min)
    /// where E_min is the minimum tabulated incoming neutrino energy, E_in is
    /// the incoming neutrino energy for the table entry, and E_out is the 
    /// out-going lepton energy. 
    /// The values in energies and zs must be sorted in ascending order (and the
    /// other tables must have been made to correspond to this). energies and zs
    /// must each contain at least three distinct values. 
    /// The s_{CC|NC}_{nu|nubar} tables must be one-dimensional datasets whose 
    /// values are the base-ten logarithms of the total cross section, in square
    /// centimeters, for the indicated current and impinging particle type. The
    /// cross section values must be tabulated for the same energies that are 
    /// contained in the energies dataset. 
    /// The dsdy_{CC|NC}_{nu|nubar} tables must be two-dimensional datasets
    /// whose values are the base-ten logarithms of the differential cross 
    /// section with respect to y, in square centimeters, for the indicated 
    /// current and impinging particle type. The dimensions of these datasets 
    /// are energy and 'z', corresponding to the values in the energies and zs
    /// datasets. Since 'z' is not well-defined for the lowest incoming neutrino 
    /// energy in the tables, E_in = E_min, the row corresponding to this energy 
    /// should contain a small value (but not negative infinity), which will 
    /// never be used directly. -50, corresponding to a cross section of 1e-50 
    /// cm^2 can be a good choice in practice. 
    void ReadHDF(const std::string& path);
    
    /// Write out an HDF5 file in the format expected by ReadHDF. 
    ///\param path the path to which the output should be written
    void WriteHDF(const std::string& path, unsigned int compressionLevel=0) const;
    
    /// \brief Returns the minimum energy in [eV]
    double GetEmin() const {return Emin;}
    /// \brief Returns the maximum energy in [eV]
    double GetEmax() const {return Emax;}
};

/// \class NeutrinoGRCrossSection
/// \brief Implements electron-antineutrino/electron scattering
class GlashowResonanceCrossSection : public NeutrinoCrossSections {
public:
  virtual ~GlashowResonanceCrossSection();

  GlashowResonanceCrossSection();
  
  GlashowResonanceCrossSection(const GlashowResonanceCrossSection&);
  GlashowResonanceCrossSection(GlashowResonanceCrossSection&&);

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
  const double fermi_scale;
  ///W mass
  const double M_W;
  ///W full decay width
  const double W_total;
  
  ///W^+ -> e^+ + \nu_e branching ratio
  static double B_Electron;
  ///W^+ -> \mu^+ + \nu_\mu branching ratio
  static double B_Muon;
  ///W^+ -> \tau^+ + \nu_\tau branching ratio
  static double B_Tau;
  ///W^+ -> hadrons branching ratio
  static double B_Hadronic;
};
	
///\brief A type alias for PDG MC numbering codes.
///The largest magnitude codes in the PDG MC numbering scheme at this time are nuclear codes
///which are "10-digit numbers ±10LZZZAAAI". Given that sign must be included, this fits within 32 
///bits with a factor of ~2 to spare.
///It is not necessary that every possible particle type be enumerated here, as non-enumerated 
///values within the range of the type may still be used. We therefore include only the most common
///codes explictly.
enum PDGCode : int32_t{
    electron=11,
    nu_e=12,
    muon=13,
    nu_mu=14,
    tau=15,
    nu_tau=16,
    ///A pseduoparticle with the average properties of a proton and neutron
    isoscalar_nucleon=81,
    proton=2212,
    neutron=2112,
    // Most nuclear targets go as 100ZZZAAA0, assuming ground state and no s quarks
    deuteron=1000010020,  // Z=1, A=2
    carbon=1000060120,    // Z=6, A=12
    oxygen=1000080160,    // Z=8, A=16
    sodium=1000110230,    // Z=11, A=23
    magnesium=1000120240, // Z=12, A=24
    aluminum=1000130270,  // Z=13, A=27
    silicon=1000140280,   // Z=14, A=28
    sulfur=1000160320,    // Z=16, A=32
    calcium=100020400,    // Z=20, A=40
    iron=1000260560,      // Z=26, A=56
    nickel=1000280580,    // Z=28, A=58
    tungsten=1000741840,  // Z=74, A=184
    lead=1000822080,      // Z=82, A=208
};
  
namespace detail{
struct PDGCodeHash{
private:
    using underlying_type=typename std::underlying_type<PDGCode>::type;
    using hash_type=std::hash<underlying_type>;
public:
    using result_type=typename hash_type::result_type;
    using argument_type=PDGCode;
    result_type operator()(argument_type arg) const{
        return hash_(arg);
    }
private:
    hash_type hash_;
};
}
  
class CrossSectionLibrary{
public:
    ///The type used to map target type codes to cross section objects
    using MapType=std::unordered_map<PDGCode,std::shared_ptr<const NeutrinoCrossSections>,detail::PDGCodeHash>;
    
    CrossSectionLibrary(){}
    CrossSectionLibrary(const MapType& crosssections):data(crosssections){}
    ///\return the requested cross section, or a null pointer if it is not found
    std::shared_ptr<const NeutrinoCrossSections> crossSectionForTarget(PDGCode target) const;
    bool hasTarget(PDGCode target) const;
    
    template <typename CrossSection>
    void addTarget(PDGCode target, CrossSection&& xs){
        if(hasTarget(target))
      throw std::runtime_error("Attempt to redefine existing target "+std::to_string(target));
        data.emplace(target, std::make_shared<CrossSection>(std::move(xs)));
    }
    void addTarget(PDGCode target, std::shared_ptr<NeutrinoCrossSections> xs){
        if(hasTarget(target))
            throw std::runtime_error("Attempt to redefine existing target "+std::to_string(target));
        data.emplace(target, xs);
    }

    int numberOfTargets();
    std::vector<PDGCode> targets();
private:
    MapType data;
};
    
CrossSectionLibrary loadDefaultCrossSections();

} // close namespace

#endif
