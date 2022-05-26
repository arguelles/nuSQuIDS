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


#ifndef NUSQUIDS_BODY_H
#define NUSQUIDS_BODY_H

#if __cplusplus < 201103L
#error C++11 compiler required. Update your compiler and use the flag -std=c++11
#endif

#include <cassert>
#include <functional>
#include <string>
#include <memory>
#include <vector>
#include <H5Ipublic.h>
#include "nuSQuIDS/marray.h"
#include "nuSQuIDS/tools.h"

namespace nusquids{

/// \class Body
/// \brief Abstract body class.
/// \details This abstract class serves as a prototype
/// of neutrino propagation environments. When implementing
/// the subclasses one must define the density() and ye()
class Body{
  protected:
    /// \brief Body parameters.
    std::vector<double> BodyParams;
    /// \brief bool that signal if the body is constant density or not
    bool is_constant_density = false;
  public:
    /// \brief Body Constructor;
    Body(){}
    virtual ~Body(){}

    /// \brief Serialization function
    virtual void Serialize(hid_t group) const=0;
    /// \brief Deserialization function
    static std::shared_ptr<Body> Deserialize(hid_t group);

    /// \class Track
    /// \brief Trajectory subclass. Specifies the trajectory
    /// to take within the body.
    class Track{
        friend class Body;
        protected:
          /// \brief Current position.
          double x;
          /// \brief Initial position.
          double xini;
          /// \brief Final position.
          double xend;
        public:
          /// \brief Default constructor.
          Track(double x,double xini, double xend): x(x), xini(xini),xend(xend) {}
          /// \brief Default constructor.
          Track(double xini, double xend): Track(xini,xini,xend) {}
          /// \brief Serialization function
          virtual void Serialize(hid_t group) const=0;
          /// \brief Deserialization function
          static std::shared_ptr<Body::Track> Deserialize(hid_t group);
          /// \brief Get track object name
          static std::string GetName() {return "BodyTrack";};
          /// \brief Sets the current position along the trajectory.
          void SetX(double y){ x = y; }
          /// \brief Returns the current position along the trajectory.
          double GetX() const { return x; }
          /// \brief Returns the initial position measured along the trajectory.
          double GetInitialX() const { return xini; }
          /// \brief Returns the final position measured along the trajectory.
          double GetFinalX() const { return xend; }
          /// \brief Reverses the track initial and final position
          void ReverseTrack() {
            double xtmp = xend;
            xend = xini;
            xini = xtmp;
          }
          /// \brief Returns parameters that define the trajectory.
          std::vector<double> GetTrackParams() const {
            std::vector<double> TrackParams{xini,xend};
            FillDerivedParams(TrackParams);
            return TrackParams;
          }
          /// Should be implemented by derived classes to append their
          /// additional parameters to TrackParams
          virtual void FillDerivedParams(std::vector<double>& TrackParams) const{}
    };
    /// \brief Return the density at a given trajectory object.
    virtual double density(const Track&) const {return 0.0;}
    /// \brief Return the electron fraction at a given trajectory object.
    virtual double ye(const Track&) const {return 1.0;}
    /// \brief Returns parameters that define the body.
    const std::vector<double>& GetBodyParams() const { return BodyParams;}
    /// \brief Returns the body identifier.
    static unsigned int GetId() {return 0;}
    /// \brief Returns the name of the body.
    static std::string GetName() {return "Body";}
    /// \brief Return true if the body is a constant density.
    virtual bool IsConstantDensity() const {return is_constant_density;}
    /// \brief Return true if the body is a constant density.
    virtual void SetIsConstantDensity(bool icd) {is_constant_density = icd;}
};

// type defining
typedef Body::Track GenericTrack;

/// \class Vacuum
/// \brief Vacuum body.
/// \details A zero density body.
class Vacuum: public Body {
  public:
    /// \brief Constructor.
    Vacuum():Body(){}

    /// \brief Serialization function
    void Serialize(hid_t group) const;
    /// \brief Deserialization function
    static std::shared_ptr<Vacuum> Deserialize(hid_t group);

    /// \class Track
    /// \brief Vacuum trajectory
    class Track: public Body::Track {
      public :
        /// \brief Rectilinear trajectory from xini to xend.
        /// @param x current position [eV^-1].
        /// @param xini Initial position [eV^-1].
        /// @param xend Final position [eV^-1].
        Track(double x,double xini,double xend):Body::Track(x,xini,xend){}
        /// \brief Rectilinear trajectory from xini to xend.
        /// @param xini Initial position [eV^-1].
        /// @param xend Final position [eV^-1].
        Track(double xini,double xend):Track(xini,xini,xend){}
        /// \brief Construct a trajectory of distance xend.
        /// @param xend Final position in eV^-1.
        /// \details In this case initial position is assumed 0.
        Track(double xend):Track(0.0,xend){}
        /// \brief Serialization function
        void Serialize(hid_t group) const;
        /// \brief Deserialization function
        static std::shared_ptr<Vacuum::Track> Deserialize(hid_t group);
        /// \brief Get track object name
        static std::string GetName() {return "VacuumTrack";}
    };
    /// \brief Returns the body identifier.
    static unsigned int GetId() {return 1;}
    /// \brief Returns the name of the body.
    static std::string GetName() {return "Vacuum";}
    /// \brief Returns the density in g/cm^3
    double density(const GenericTrack&) const;
    /// \brief Returns the electron fraction
    double ye(const GenericTrack&) const;
    /// \brief Returns true as this body has constant density by definition
    bool IsConstantDensity() const;
};

/// \class ConstantDensity
/// \brief Constant density body.
/// \details A body with constant density
class ConstantDensity: public Body{
  private:
    /// \brief Constant density in g/cm^3
    const double constant_density;
    /// \brief Constant electron fraction
    const double constant_ye;
  public:
    /// \brief Constructor
    /// @param density Density in g/cm^3.
    /// @param ye electron fraction.
    ConstantDensity(double density,double ye);

    /// \brief Serialization function
    void Serialize(hid_t group) const;
    /// \brief Deserialization function
    static std::shared_ptr<ConstantDensity> Deserialize(hid_t group);

    /// \class Track
    /// \brief Constant density trajectory
    class Track: public Body::Track{
      public :
        /// \brief Rectilinear trajectory from xini to xend.
        /// @param x current position [eV^-1].
        /// @param xini Initial position [eV^-1].
        /// @param xend Final position [eV^-1].
        Track(double x,double xini,double xend):Body::Track(x,xini,xend){}
        /// \brief Construct a trajectory between an initial position and final position.
        /// @param xini Initial position in eV^-1.
        /// @param xend Final position in eV^-1.
        Track(double xini,double xend):Track(xini,xini,xend){}
        /// \brief Construct a trajectory of distance xend.
        /// @param xend Final position in eV^-1.
        /// \details In this case initial position is assumed 0.
        Track(double xend):Track(0.0,xend){}
        /// \brief Serialization function
        void Serialize(hid_t group) const;
        /// \brief Deserialization function
        static std::shared_ptr<ConstantDensity::Track> Deserialize(hid_t group);
        /// \brief Get track object name
        static std::string GetName() {return "ConstantDensityTrack";}
    };

    /// \brief Returns the body identifier.
    static unsigned int GetId() {return 2;}
    /// \brief Returns the name of the body.
    static std::string GetName() {return "ConstantDensity";}
    /// \brief Returns the density in g/cm^3
    double density(const GenericTrack&) const;
    /// \brief Returns the electron fraction
    double ye(const GenericTrack&) const;
    /// \brief Returns true as this body has constant density by definition
    bool IsConstantDensity() const;
};

/// \class VariableDensity
/// \brief Variable user defined density body.
/// \details The user provides arrays which contain the body
/// information.
class VariableDensity: public Body{
  private:
    /// \brief Array of length \c arraysize containing spline position nodes.
    std::vector<double> x_arr;
    /// \brief Array of length \c arraysize containing the density at each position nodes.
    std::vector<double> density_arr;
    /// \brief Array of length \c arraysize containing the electron fraction at each position nodes.
    std::vector<double> ye_arr;
    /// \brief Size of array used to create spline.
    unsigned int arraysize;

    /// \brief Minimum value of \c x_arr.
    double x_max;
    /// \brief Maximum value of \c x_arr.
    double x_min;

    /// \brief Density spline
    AkimaSpline inter_density;
    /// \brief Electron fraction spline
    AkimaSpline inter_ye;
  public:
    /// \brief Constructor.
    /// @param x Vector containing position nodes in cm.
    /// @param density Density, in g/cm^3, at each of the nodes.
    /// @param ye Electron fraction at each of the nodes.
    /// \pre All input vectors must be of equal size.
    VariableDensity(std::vector<double> x,std::vector<double> density,std::vector<double> ye);
    /// \brief Destructor.
    ~VariableDensity();

    /// \brief Serialization function
    void Serialize(hid_t group) const;
    /// \brief Deserialization function
    static std::shared_ptr<VariableDensity> Deserialize(hid_t group);

    /// \class Track
    /// \brief Variable density trajectory
    class Track: public Body::Track{
      public :
        /// \brief Rectilinear trajectory from xini to xend.
        /// @param x current position [eV^-1].
        /// @param xini Initial position [eV^-1].
        /// @param xend Final position [eV^-1].
        Track(double x,double xini,double xend):Body::Track(x,xini,xend){}
        /// \brief Construct a trajectory between an initial position and final position.
        /// @param xini Initial position in eV^-1.
        /// @param xend Final position in eV^-1.
        Track(double xini,double xend):Track(xini,xini,xend){}
        /// \brief Construct a trajectory of distance xend.
        /// @param xend Final position in eV^-1.
        /// \details In this case initial position is assumed 0.
        Track(double xend):Track(0.0,xend){}
        /// \brief Serialization function
        void Serialize(hid_t group) const;
        /// \brief Deserialization function
        static std::shared_ptr<VariableDensity::Track> Deserialize(hid_t group);
        /// \brief Get track object name
        static std::string GetName() {return "VariableDensityTrack";}

    };
    /// \brief Returns the body identifier.
    static unsigned int GetId() {return 3;}
    /// \brief Returns the name of the body.
    static std::string GetName() {return "VariableDensity";}
    /// \brief Returns the density in g/cm^3
    double density(const GenericTrack&) const;
    /// \brief Returns the electron fraction
    double ye(const GenericTrack&) const;
};

/// \class Earth
/// \brief Earth model based on PREM.
class Earth: public Body{
  private:
    /// \brief Radius of the Earth.
    double radius;
    /// \brief Earth radius position array
    std::vector<double> earth_radius;
    /// \brief Earth density array
    std::vector<double> earth_density;
    /// \brief Earth electron fraction array
    std::vector<double> earth_ye;
    /// \brief Data arrays size
    unsigned int arraysize;

    /// \brief Density spline
    AkimaSpline inter_density;
    /// \brief Electron fraction spline
    AkimaSpline inter_ye;

    /// \brief Minimum radius.
    double x_radius_min;
    /// \brief Maximum radius.
    double x_radius_max;
    /// \brief Density at minimum radius.
    double x_rho_min;
    /// \brief Density at maximum radius.
    double x_rho_max;
    /// \brief Electron fraction at minimum radius.
    double x_ye_min;
    /// \brief Electron fraction at maximum radius.
    double x_ye_max;
  public:
    /// \brief Default constructor using supplied PREM.
    Earth();
    /// \brief Constructor from a user supplied Earth model.
    /// @param earthmodel Path to the Earth model file.
    /// \details The input file should have three columns.
    /// The first one must run from zero to one representing
    /// the center and surface of the Earth respectively. The
    /// second column must contain the Earth density in g/cm^3 at
    /// a given position, while the third column must contain
    /// the electron fraction.
    Earth(std::string earthmodel);
    /// \brief Constructor in which the user provides, as vectors, the
    /// Earth properties.
    /// @param x Vector containing position nodes in cm.
    /// @param rho Density, in g/cm^3, at each of the nodes.
    /// @param ye Electron fraction at each of the nodes.
    /// \pre All input vectors must be of equal size.
    Earth(std::vector<double> x,std::vector<double> rho,std::vector<double> ye);
    /// \brief Destructor.
    ~Earth();

    /// \brief Serialization function
    void Serialize(hid_t group) const;
    /// \brief Deserialization function
    static std::shared_ptr<Earth> Deserialize(hid_t group);

    /// \class Track
    /// \brief Earth trajectory
    class Track: public Body::Track{
      private:
        /// \brief Baseline of the neutrino experiment in natural units.
        const double baseline;
      public :
        /// \brief Rectilinear trajectory from xini to xend.
        /// @param x current position [eV^-1].
        /// @param xini Initial position [eV^-1].
        /// @param xend Final position [eV^-1].
        Track(double x,double xini,double xend,double baseline):Body::Track(x,xini,xend),baseline(baseline){}
        /// \brief Construct a trajectory between an initial position and final position.
        /// @param xini Initial position in eV^-1.
        /// @param xend Final position in eV^-1.
        /// @param baseline Baseline of experiment in eV^-1.
        Track(double xini,double xend,double baseline):Track(xini,xini,xend,baseline){}
        /// \brief Construct a trajectory where the neutrino travels a distance given by baseline.
        /// @param baseline Traverse distance in eV^-1.
        /// \details In this case \c xini = 0, \c xend = \c baseline.
        Track(double baseline):Track(0.,baseline,baseline){}
        /// \brief Returns the neutrino baseline in natural units.
        double GetBaseline() const {return baseline;}
        virtual void FillDerivedParams(std::vector<double>& TrackParams) const;
        /// \brief Serialization function
        void Serialize(hid_t group) const;
        /// \brief Deserialization function
        static std::shared_ptr<Earth::Track> Deserialize(hid_t group);
        /// \brief Get track object name
        static std::string GetName() {return "EarthTrack";}
    };
    /// \brief Returns the body identifier.
    static unsigned int GetId() {return 4;}

    /// \brief Returns the name of the body.
    static std::string GetName() {return "Earth";}

    /// \brief Returns the density in g/cm^3
    double density(const GenericTrack&) const;
    /// \brief Returns the electron fraction
    double ye(const GenericTrack&) const;

    /// \brief Returns the radius of the Earth in natural units.
    double GetRadius() const {return radius;}
};

/// \class Sun
/// \brief A model of the Sun.
class Sun: public Body{
  private:
    /// \brief Array that contains all the Solar model parameters.
    marray<double,2> sun_model;
    /// \brief Array that contains all the Solar electron density parameters.
    marray<double,2> sun_model_nele;
    /// \brief Array of length \c arraysize containing spline position nodes.
    std::vector<double> sun_radius;
    /// \brief Array of length \c arraysize containing the density at position nodes.
    std::vector<double> sun_density;
    /// \brief Array of length \c arraysize containing the hydrogen fraction at position nodes.
    std::vector<double> sun_xh;

    /// \brief Size of \c sun_radius array.
    unsigned int arraysize;

    /// \brief Radius of the Sun.
    double radius;

    /// \brief Density spline
    AkimaSpline inter_density;
    /// \brief Hydrogen fraction spline
    AkimaSpline inter_xh;

    /// \brief Returns the density in g/cm^3 at a given radius fraction x
    /// @param x Radius fraction: 0:center, 1:surface.
    double rdensity(double x) const;
    /// \brief Returns the electron fraction at a given radius fraction x
    /// @param x Radius fraction: 0:center, 1:surface.
    double rxh(double x) const;
  public:
    /// \brief Detault constructor.
    Sun();
    /// \brief Constructor in which the user provides, as vectors, the
    /// Sun properties.
    /// @param x Vector containing position nodes in cm.
    /// @param rho Density, in g/cm^3, at each of the nodes.
    /// @param xh Hydrogen fraction at each point.
    /// \pre All input vectors must be of equal size.
    Sun(std::vector<double> x,std::vector<double> rho,std::vector<double> xh);
    /// \brief Constructor from a user supplied solar model.
    /// @param sunmodel Path to the Sun model file.
    /// \details The input file should have the same columns as John Bahcalls solar model file.
    /// The second column one must run from zero to one representing
    /// the center and surface of the Sun respectively. The
    /// fourth column must contain the Sun density in g/cm^3 at
    /// a given position, while the seven column must contain
    /// the hydrogen fraction which is related to the electron fraction by ye = 0.5*(1.0+rxh(r)).
    /// Internally we assume that the solar radius is 695980.0 kilometers.
    Sun(std::string sunmodel);


    /// \brief Destructor
    ~Sun();

    /// \brief Serialization function
    void Serialize(hid_t group) const;
    /// \brief Deserialization function
    static std::shared_ptr<Sun> Deserialize(hid_t group);

    /// \class Track
    /// \brief Sun trajectory
    class Track: public Body::Track{
      public :
        /// \brief Rectilinear trajectory from xini to xend.
        /// @param x current position [eV^-1].
        /// @param xini Initial position [eV^-1].
        /// @param xend Final position [eV^-1].
        Track(double x,double xini,double xend):Body::Track(x,xini,xend){}
        /// \brief Construct a trajectory between an initial position and final position.
        /// @param xini Initial position in eV^-1.
        /// @param xend Final position in eV^-1.
        /// \details The trajectory is measured from the sun center which is set to zero.
        Track(double xini,double xend):Track(xini,xini,xend){}
        /// \brief Construct a trajectory from the sun center to a final position.
        /// @param xend Final position in eV^-1.
        /// \details The trajectory is measured from the sun center which is set to zero.
        Track(double xend):Track(0.,xend){}
        /// \brief Serialization function
        void Serialize(hid_t group) const;
        /// \brief Deserialization function
        static std::shared_ptr<Sun::Track> Deserialize(hid_t group);
        /// \brief Get track object name
        static std::string GetName() {return "SunTrack";}
    };

    /// \brief Returns the body identifier.
    static unsigned int GetId() {return 5;}

    /// \brief Returns the name of the body.
    static std::string GetName() {return "Sun";}

    /// \brief Returns the density in g/cm^3
    double density(const GenericTrack&) const;
    /// \brief Returns the electron fraction
    double ye(const GenericTrack&) const;

    /// \brief Returns the radius of the Sun in natural units.
    double GetRadius() const {return radius;}
};

/// \class SunASnu
/// \brief A model of the Sun with atmospheric solar neutrinos geometry.
class SunASnu: public Body{
  private:
    /// \brief Array that contains all the Solar model parameters.
    marray<double,2> sun_model;
    /// \brief Array of length \c arraysize containing spline position nodes.
    std::vector<double> sun_radius;
    /// \brief Array of length \c arraysize containing the density at position nodes.
    std::vector<double> sun_density;
    /// \brief Array of length \c arraysize containing the hydrogen fraction at position nodes.
    std::vector<double> sun_xh;

    /// \brief Size of \c sun_radius array.
    unsigned int arraysize;

    /// \brief Radius of the Sun.
    double radius;

    /// \brief Density spline
    AkimaSpline inter_density;
    /// \brief Hydrogen fraction spline
    AkimaSpline inter_xh;

    /// \brief Returns the density in g/cm^3 at a given radius fraction x
    /// @param x Radius fraction: 0:center, 1:surface.
    double rdensity(double x) const;
    /// \brief Returns the electron fraction at a given radius fraction x
    /// @param x Radius fraction: 0:center, 1:surface.
    double rxh(double x) const;
  public:
    /// \brief Detault constructor.
    SunASnu();
    /// \brief Constructor from a user supplied solar model.
    /// @param sunmodel Path to the Sun model file.
    /// \details The input file should have the same columns as John Bahcalls solar model file.
    /// The second column one must run from zero to one representing
    /// the center and surface of the Sun respectively. The
    /// fourth column must contain the Sun density in g/cm^3 at
    /// a given position, while the seven column must contain
    /// the hydrogen fraction which is related to the electron fraction by ye = 0.5*(1.0+rxh(r)).
    /// Internally we assume that the solar radius is 695980.0 kilometers.
    SunASnu(std::string sunmodel);
    /// \brief Constructor in which the user provides, as vectors, the
    /// Sun properties.
    /// @param x Vector containing position nodes in cm.
    /// @param rho Density, in g/cm^3, at each of the nodes.
    /// @param xh Hydrogen fraction at each point.
    /// \pre All input vectors must be of equal size.
    SunASnu(std::vector<double> x,std::vector<double> rho,std::vector<double> xh);

    ~SunASnu();

    /// \brief Serialization function
    void Serialize(hid_t group) const;
    /// \brief Deserialization function
    static std::shared_ptr<SunASnu> Deserialize(hid_t group);

    /// \class Track
    /// \brief SunASnu trajectory
    class Track: public Body::Track{
      friend class SunASnu;
      private:
        /// \brief Radius of the Sun.
        double radius_nu;
        /// \brief Impact parameter.
        double b_impact;
      public:
        /// \brief Rectilinear trajectory from xini to xend.
        /// @param x current position [eV^-1].
        /// @param xini Initial position [eV^-1].
        /// @param b_impact impact parameter in eV^-1.
        Track(double x,double xini,double b_impact);
        /// \brief Construct a trajectory between an initial position and final position.
        /// @param xini Initial position in eV^-1.
        /// @param b_impact impact parameter in eV^-1.
        /// \details The trajectory baseline is determined by the impact parameter and starts
        /// at \c xini, and ends when the neutrino exits the sun.
        Track(double xini,double b_impact):Track(xini,xini,b_impact){}
        /// \brief Construct a for a given impact parameter.
        /// @param b_impact_ impact parameter in eV^-1.
        /// \details The trajectory baseline is determined by the impact parameter and starts
        /// at \c xini = 0, and ends when the neutrino exits the sun.
        Track(double b_impact_):Track(0.0,b_impact_){}
        virtual void FillDerivedParams(std::vector<double>& TrackParams) const;
        /// \brief Serialization function
        void Serialize(hid_t group) const;
        /// \brief Deserialization function
        static std::shared_ptr<SunASnu::Track> Deserialize(hid_t group);
        /// \brief Get track object name
        static std::string GetName() {return "SunASnuTrack";}
    };

    /// \brief Returns the body identifier.
    static unsigned int GetId() {return 6;}

    /// \brief Returns the name of the body.
    static std::string GetName() {return "SunASnu";}

    /// \brief Returns the density in g/cm^3
    double density(const GenericTrack&) const;
    /// \brief Returns the electron fraction
    double ye(const GenericTrack&) const;

    /// \brief Returns the radius of the Sun in natural units.
    double GetRadius() const {return radius;}
};

/// \class EarthAtm
/// \brief A model of the Earth with atmospheric neutrinos geometry.
class EarthAtm: public Body{
  protected:
    /// \brief Radius of the Earth.
    double radius;
    /// \brief Height of the atmosphere.
    double atm_height;
    /// \brief Radius of the Earth plus atmosphere.
    double earth_with_atm_radius;

    /// \brief Earth radius position array
    std::vector<double> earth_radius;
    /// \brief Earth density array
    std::vector<double> earth_density;
    /// \brief Earth electron fraction array
    std::vector<double> earth_ye;
    /// \brief Data arrays size
    unsigned int arraysize;
  
    /// \brief Density spline
	AkimaSpline inter_density;
    /// \brief Electron fraction spline
    AkimaSpline inter_ye;

    /// \brief Minimum radius.
    double x_radius_min;
    /// \brief Maximum radius.
    double x_radius_max;
    /// \brief Density at minimum radius.
    double x_rho_min;
    /// \brief Density at maximum radius.
    double x_rho_max;
    /// \brief Electron fraction at minimum radius.
    double x_ye_min;
    /// \brief Electron fraction at maximum radius.
    double x_ye_max;
  public:
    /// \brief Default constructor using supplied PREM.
    EarthAtm();
    /// \brief Constructor from a user supplied Earth model.
    /// @param earthmodel Path to the Earth model file.
    /// \details The input file should have three columns.
    /// The first one must run from zero to one representing
    /// the center and surface of the Earth respectively. The
    /// second column must contain the Earth density in g/cm^3 at
    /// a given position, while the third column must contain
    /// the electron fraction.
    EarthAtm(std::string earthmodel);
    /// \brief Constructor in which the user provides, as vectors, the
    /// Earth properties.
    /// @param x Vector containing position nodes in cm.
    /// @param rho Density, in g/cm^3, at each of the nodes.
    /// @param ye Electron fraction at each of the nodes.
    /// \pre All input vectors must be of equal size.
    EarthAtm(std::vector<double> x,std::vector<double> rho,std::vector<double> ye);
    ~EarthAtm();

    /// \brief Serialization function
    void Serialize(hid_t group) const;
    /// \brief Deserialization function
    static std::shared_ptr<EarthAtm> Deserialize(hid_t group);

    /// \class Track
    /// \brief EarthAtm trajectory
    class Track: public Body::Track{
      friend class EarthAtm;
      private:
        /// \brief Cosine of the zenith angle.
        double cosphi;
        /// \brief Radius of the Earth.
        double radius_nu;
        /// \brief Height of the atmosphere.
        double atmheight;
        /// \brief Baseline.
        double L;

        ///Private constructor which does not initialize all members.
        ///Intended for use only by makeWithCosine.
        Track();
      public :
        /// \brief Construct trajectory.
        /// @param x_ current position [eV^-1].
        /// @param phi Zenith angle in radians.
        Track(double x_,double phi):Track(phi){x=x_; assert(x >= 0);};
        /// \brief Construct trajectory.
        /// @param phi Zenith angle in radians.
        Track(double phi);
        /// \brief Returns the neutrino baseline in natural units.
        double GetBaseline() const {return L;}
        virtual void FillDerivedParams(std::vector<double>& TrackParams) const;
        ///Construct a track with the cosine of the zenith angle
        static Track makeWithCosine(double cosphi);
        /// \brief Serialization function
        void Serialize(hid_t group) const;
        /// \brief Deserialization function
        static std::shared_ptr<EarthAtm::Track> Deserialize(hid_t group);
        /// \brief Get track object name
        static std::string GetName() {return "EarthAtmTrack";}
    };
    /// \brief Returns the body identifier.
    static unsigned int GetId() {return 7;}
    /// \brief Returns the name of the body.
    static std::string GetName() {return "EarthAtm";}
    /// \brief Returns the density in g/cm^3
    double density(const GenericTrack&) const;
    /// \brief Returns the electron fraction
    double ye(const GenericTrack&) const;
    /// \brief Returns the radius of the Earth in km.
    double GetRadius() const {return radius;}
    /// \brief Returns the neutrino production altitude in km.
    double GetProductionAltitude() const {return atm_height;}
};

  // type defining
  typedef Body::Track Track;

  // registration body
  namespace detail{
    class registerBody{
      public:
        registerBody(const std::string& name,std::function<std::shared_ptr<Body>(hid_t)> fdeserialize);
    };
    class registerTrack{
      public:
        registerTrack(const std::string& name,std::function<std::shared_ptr<Track>(hid_t)> fdeserialize);
    };
  }

  std::function<std::shared_ptr<Body>(hid_t)> GetBodyDeserializer(std::string);
  std::function<std::shared_ptr<Track>(hid_t)> GetTrackDeserializer(std::string);


// body registration
#define ASW(a, b) a ## b

#define NUSQUIDS_REGISTER_BODY(classname) \
namespace{ nusquids::detail::registerBody ASW(body_registerer,classname)(classname::GetName(),classname::Deserialize); nusquids::detail::registerTrack ASW(track_registerer,classname)(classname::Track::GetName(),classname::Track::Deserialize);}

} // close namespace

#endif
