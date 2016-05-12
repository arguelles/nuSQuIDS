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

/* This is a simple example where we modified the earth object, constructing a new object
 * called EarthMod, this very simple modification consist in including 3 parameters that weight 
 * density in every layer of the standard earth model.
 */


namespace nusquids{


class EarthMod: public Body{
private:
  /// \brief Radius of the EarthMod.
  double radius;
    
  /// \brief Density gsl spline.
  gsl_spline * inter_density;
  /// \brief Density gsl spline auxiliary pointer.
  gsl_interp_accel * inter_density_accel;
  /// \brief Electron fraction gsl spline.
  gsl_spline * inter_ye;
  /// \brief Electron fraction gsl spline auxiliary pointer.
  gsl_interp_accel * inter_ye_accel;
  double earth_with_atm_radius;
  double atm_height;


  marray<double,2> earth_model; 

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
  EarthMod():
  EarthMod(EARTH_MODEL_LOCATION,1,1,1){}
  /// \brief Constructor from a user supplied EarthMod model.
  /// @param earthmodel Path to the EarthMod model file.
  /// \details The input file should have three columns.
  /// The first one must run from zero to one representing
  /// the center and surface of the EarthMod respectively. The
  /// second column must contain the EarthMod density in gr/cm^3 at
  /// a given position, while the third column must contain
  /// the electron fraction.
  EarthMod(std::string earthmodel, double, double , double);
  //EarthMod(double, double , double);
  /// \brief Destructor.
  ~EarthMod();
  //This function sets the values that weight the different layers in the of the PREM model
  void  Mod(double frho1, double frho2, double frho3);  
  /// \class Track
  /// \brief EarthMod trajectory
  class Track: public Body::Track{
    friend class EarthMod;
  private:
    /// \brief Zenith angle in radians.
    double phi;
    /// \brief Cosine of the zenith angle.
    double cosphi;
    /// \brief Radius of the Earth.
    double radius_nu;
    /// \brief Height of the atmosphere.
    double atmheight;
    /// \brief Baseline.
    double L;
  public :
    /// \brief Construct trajectory.
    /// @param phi Zenith angle in radians.
    Track(double phi);
    /// \brief Returns the neutrino baseline in natural units.
    double GetBaseline() const {return L;}
  };
  
  /// \brief Returns the density in g/cm^3
  double density(const GenericTrack&) const;
  /// \brief Returns the electron fraction
  double ye(const GenericTrack&) const;
  
  /// \brief Returns the radius of the EarthMod in natural units.
  double GetRadius() const {return radius;}
};
  

  

}



#endif
