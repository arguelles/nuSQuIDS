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


#include "body.h"

// Macros
#define SQR(x)      ((x)*(x))                        // x^2

namespace nusquids{

static squids::Const param;

/*
----------------------------------------------------------------------
         VACUUM CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// track constructor
Vacuum::Track::Track(double xini, double xend):Body::Track(xini,xend)
        {
            x = xini;
            TrackParams = {xini,xend};
        }

double Vacuum::density(const GenericTrack& track_input) const{
            return 0.0;
        }

double Vacuum::ye(const GenericTrack& track_input) const{
            return 1.0;
        }

/*
----------------------------------------------------------------------
         ConstantDensity CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
ConstantDensity::ConstantDensity(double constant_density,double constant_ye):Body(2,"ConstantDensity"),
                                                                             constant_density(constant_density),
                                                                             constant_ye(constant_ye)
        {
            BodyParams = {constant_density, constant_ye};
        }

// track constructor
ConstantDensity::Track::Track(double xini, double xend):Body::Track(xini,xend)
        {
            x = xini;
            TrackParams = {xini,xend};
        }

double ConstantDensity::density(const GenericTrack& track_input) const
        {
            return constant_density;
        }

double ConstantDensity::ye(const GenericTrack& track_input) const
        {
            return constant_ye;
        }

/*
----------------------------------------------------------------------
         VariableDensity CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
VariableDensity::VariableDensity(std::vector<double> x_input,std::vector<double> density_input,std::vector<double> ye_input):Body(3,"VariableDensity")
        {
            assert("nuSQUIDS::Error::VariableDensityConstructor: Invalid array sizes." && x_input.size() == density_input.size() && x_input.size() == ye_input.size());
            arraysize = x_input.size();

            x_min = x_input.front();
            x_max = x_input.back();

            x_arr = new double[arraysize];
            density_arr = new double[arraysize];
            ye_arr = new double [arraysize];

            for(unsigned int i = 0; i < arraysize; i++){
              x_arr[i] = x_input[i];
              density_arr[i] = density_input[i];
              ye_arr[i] = ye_input[i];
            }

            inter_density = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_density_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_density,x_arr,density_arr,arraysize);

            inter_ye = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_ye_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_ye,x_arr,ye_arr,arraysize);

            for(double xx : x_input)
              BodyParams.push_back(xx);
            for(double rho : density_input)
              BodyParams.push_back(rho);
            for(double ye : ye_input)
              BodyParams.push_back(ye);
        }

// track constructor
VariableDensity::Track::Track(double xini, double xend):Body::Track(xini,xend)
        {
            x = xini;
            TrackParams = {xini,xend};
        }

double VariableDensity::density(const GenericTrack& track_input) const
        {
          double x = track_input.GetX()/param.cm;
          if (x < x_min or x > x_max ){
              return 0;
          } else {
              return gsl_spline_eval(inter_density,x,inter_density_accel);
          }
        }
double VariableDensity::ye(const GenericTrack& track_input) const
        {
          double x = track_input.GetX()/param.cm;
          if (x < x_min or x > x_max ){
              return 0;
          } else {
              return gsl_spline_eval(inter_ye,x,inter_ye_accel);
          }
        }

/*
----------------------------------------------------------------------
         Earth CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
Earth::Earth():Earth(static_cast<std::string>(EARTH_MODEL_LOCATION))
        {
        }

// track constructor
Earth::Track::Track(double xini, double xend,double baseline): Body::Track(xini,xend),baseline(baseline)
        {
            x = xini;
            TrackParams = {xini,xend,baseline};
        }

double Earth::density(const GenericTrack& track_input) const
        {
            //std::shared_ptr<const Earth::Track> track_earth = std::static_pointer_cast<const Earth::Track >(track_input);
            const Earth::Track& track_earth = static_cast<const Earth::Track&>(track_input);
            double xkm = track_earth.GetX()/param.km;
            double r = sqrt(SQR(radius)+SQR(xkm)-(track_earth.GetBaseline()/param.km)*xkm);

            if ( r/radius < x_radius_min ){
              return x_rho_min;
            }
            else if ( r/radius > x_radius_max ) {
              return x_rho_max;
            }
            else {
              return gsl_spline_eval(inter_density,r/radius,inter_density_accel);
            }
        }

double Earth::ye(const GenericTrack& track_input) const
        {
            //std::shared_ptr<const Earth::Track> track_earth = std::static_pointer_cast<const Earth::Track >(track_input);
            const Earth::Track& track_earth = static_cast<const Earth::Track&>(track_input);
            double xkm = track_earth.GetX()/param.km;
            double r = sqrt(SQR(radius)+SQR(xkm)-(track_earth.GetBaseline()/param.km)*xkm);

            if ( r/radius < x_radius_min ){
              return x_ye_min;
            }
            else if ( r/radius > x_radius_max ) {
              return x_ye_max;
            }
            else {
              return gsl_spline_eval(inter_ye,r/radius,inter_ye_accel);
            }
        }

Earth::Earth(std::string filepath):Body(4,"Earth")
        {
          // The Input file should have the radius specified from 0 to 1.
          // where 0 is the center of the Earth and 1 is the surface.
            radius = 6371.0; // [km]

            marray<double,2> earth_model = quickread(filepath);
            size_t arraysize = earth_model.extent(0);

            double earth_radius[arraysize];
            double earth_density[arraysize];
            double earth_ye[arraysize];

            for (unsigned int i=0; i < arraysize;i++){
                earth_radius[i] = earth_model[i][0];
                earth_density[i] = earth_model[i][1];
                earth_ye[i] = earth_model[i][2];
            }

            x_radius_min = earth_radius[0];
            x_radius_max = earth_radius[arraysize-1];
            x_rho_min = earth_density[0];
            x_rho_max = earth_density[arraysize-1];
            x_ye_min = earth_ye[0];
            x_ye_max = earth_ye[arraysize-1];

            inter_density = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_density_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_density,earth_radius,earth_density,arraysize);

            inter_ye = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_ye_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_ye,earth_radius,earth_ye,arraysize);
        }

Earth::~Earth(){
  gsl_spline_free(inter_density);
  gsl_interp_accel_free(inter_density_accel);
  gsl_spline_free(inter_ye);
  gsl_interp_accel_free(inter_ye_accel);
}

/*
----------------------------------------------------------------------
         SUN CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
Sun::Sun():Body(5,"Sun")
        {
            radius = 695980.0*param.km;

            // import sun model
            sun_model = quickread(SUN_MODEL_LOCATION);
            arraysize = sun_model.extent(0);

            sun_radius = new double[arraysize];
            sun_density = new double[arraysize];
            sun_xh = new double[arraysize];

            for (unsigned int i=0; i < arraysize;i++){
                sun_radius[i] = sun_model[i][1];
                sun_density[i] = sun_model[i][3];
                sun_xh[i] = sun_model[i][6];
            }

            inter_density = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_density_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_density,sun_radius,sun_density,arraysize);

            inter_rxh = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_rxh_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_rxh,sun_radius,sun_xh,arraysize);
        }
// track constructor
Sun::Track::Track(double xini, double xend):Body::Track(xini,xend)
        {
            x = xini;
            TrackParams = {xini,xend};
        }

double Sun::rdensity(double x) const{
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_density[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_density,x,inter_density_accel);
            }
        }

double Sun::rxh(double x) const{
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_xh[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_rxh,x,inter_rxh_accel);
            }
        }

double Sun::density(const GenericTrack& track_input) const
        {
            double r = track_input.GetX()/(radius);
            return rdensity(r);
        }
double Sun::ye(const GenericTrack& track_input) const
        {
            double r = track_input.GetX()/(radius);
            return 0.5*(1.0+rxh(r));
        }

Sun::~Sun(){
  free(sun_radius);
  free(sun_density);
  free(sun_xh);
  //free(sun_nele_radius);
  //free(sun_nele);
  gsl_spline_free(inter_density);
  gsl_interp_accel_free(inter_density_accel);
  gsl_spline_free(inter_rxh);
  gsl_interp_accel_free(inter_rxh_accel);
  //free(inter_nele);
  //free(inter_nele_accel);
}
/*
----------------------------------------------------------------------
         SUN ASNU CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
SunASnu::SunASnu():Body(6,"SunASnu")
        {
            radius = 694439.0*param.km;

            sun_model = quickread(SUN_MODEL_LOCATION);
            arraysize = sun_model.extent(0);

            sun_radius = new double[arraysize];
            sun_density = new double[arraysize];
            sun_xh = new double[arraysize];

            for (unsigned int i=0; i < arraysize;i++){
                sun_radius[i] = sun_model[i][1];
                sun_density[i] = sun_model[i][3];
                sun_xh[i] = sun_model[i][6];
            }

            inter_density = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_density_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_density,sun_radius,sun_density,arraysize);

            inter_rxh = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_rxh_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_rxh,sun_radius,sun_xh,arraysize);
        }
// track constructor
SunASnu::Track::Track(double xini, double b_impact):
  Body::Track(xini,xini),
  radius_nu(694439.0*param.km),
  b_impact(b_impact)
        {
            x = xini;
            xend = 2.0*sqrt(SQR(radius_nu)+SQR(b_impact));
            TrackParams = {xini,xend,b_impact};
        }

double SunASnu::rdensity(double x) const{
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_density[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_density,x,inter_density_accel);
            }
        }

double SunASnu::rxh(double x) const{
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_xh[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_rxh,x,inter_rxh_accel);
            }
        }

double SunASnu::density(const GenericTrack& track_input) const
        {
            //std::shared_ptr<const SunASnu::Track> track_sunasnu = std::static_pointer_cast<const SunASnu::Track >(track_input);
            const SunASnu::Track& track_sunasnu = static_cast<const SunASnu::Track&>(track_input);
            double x = track_sunasnu.GetX();
            double b = track_sunasnu.b_impact;

            double r = sqrt(SQR(radius)+SQR(x)-2.0*x*sqrt(SQR(radius)-SQR(b)))/radius;

            return rdensity(r);
        }

double SunASnu::ye(const GenericTrack& track_input) const
        {
            //std::shared_ptr<const SunASnu::Track> track_sunasnu = std::static_pointer_cast<const SunASnu::Track >(track_input);
            const SunASnu::Track& track_sunasnu = static_cast<const SunASnu::Track&>(track_input);
            double x = track_sunasnu.GetX();
            double b = track_sunasnu.b_impact;
            double r = sqrt(SQR(radius)+SQR(x)-2.0*x*sqrt(SQR(radius)-SQR(b)))/radius;
            return 0.5*(1.0+rxh(r));
        }

SunASnu::~SunASnu(){
  free(sun_radius);
  free(sun_density);
  free(sun_xh);
  gsl_spline_free(inter_density);
  gsl_interp_accel_free(inter_density_accel);
  gsl_spline_free(inter_rxh);
  gsl_interp_accel_free(inter_rxh_accel);
}

/*
----------------------------------------------------------------------
         EARTHATM CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
EarthAtm::EarthAtm():EarthAtm(EARTH_MODEL_LOCATION)
        {
        }

// track constructor
EarthAtm::Track::Track(double phi_input):Body::Track(0,0)
        {
            radius_nu = 6371.0*param.km;
            radius_nu = 6369.0*param.km;
            atmheight = 100.0*param.km;
            atmheight = 22.;

            phi = phi_input;
            cosphi = cos(phi);

            /*
            if(cosphi<=0.0){
                L = 2.0*radius_nu*std::abs(cosphi);
            } else {
                L = atmheight/std::abs(cosphi);
            }
            */

            double R = radius_nu;
            double r = atmheight;
            double mm = tan(phi);

            if(cosphi<=0.0){
                L = sqrt(((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*r*R + 2.0*R*
                    (R + sqrt((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*
                    r*R + R*R)))/(1.0 + mm*mm));
            } else {
                L = sqrt(((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*r*R + 2.0*R*
                    (R - sqrt((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*
                    r*R + R*R)))/(1.0 + mm*mm));
            }

            x = 0.0;
            xini = 0.0;
            xend = L;

            #ifdef EarthAtm_DEBUG
                cout << "== Init Track ==" << endl;
                cout << " phi = " << phi <<
                ", cos(phi) = " << cosphi <<
                ", L = " << radius_nu/param.km << endl;
                cout << "==" << endl;
            #endif

            TrackParams = {xini,xend,phi_input};
        }

double EarthAtm::density(const GenericTrack& track_input) const
        {
            const EarthAtm::Track& track_earthatm = static_cast<const EarthAtm::Track&>(track_input);
            double xkm = track_earthatm.GetX()/param.km;

            double r = sqrt(SQR(earth_with_atm_radius) + SQR(xkm) - (track_earthatm.L/param.km)*xkm);

            #ifdef EarthAtm_DEBUG
            cout << "r : " << r << " L : " << (track_earthatm->L/param.km)
                 << " x : " << xkm << " R : " << radius << endl;
            #endif

            double rel_r = r/earth_with_atm_radius;

            if ( rel_r < x_radius_min ){
              return x_rho_min;
            }
            else if ( rel_r > x_radius_max and rel_r < radius/earth_with_atm_radius) {
              return x_rho_max;
            }
            else if ( rel_r > radius/earth_with_atm_radius ) {
                double h = atm_height*(rel_r - radius/earth_with_atm_radius);
                double h0 = 25.0;
                return 1.05*exp(-h/h0);
            } else {
              return gsl_spline_eval(inter_density,r/radius,inter_density_accel);
            }
        }

double EarthAtm::ye(const GenericTrack& track_input) const
        {
            const EarthAtm::Track& track_earthatm = static_cast<const EarthAtm::Track&>(track_input);
            double xkm = track_earthatm.GetX()/param.km;
            double r = sqrt(SQR(earth_with_atm_radius) + SQR(xkm) - (track_earthatm.L/param.km)*xkm);

            double rel_r = r/earth_with_atm_radius;
            if ( rel_r < x_radius_min ){
              return x_ye_min;
            }
            else if ( rel_r > x_radius_max and rel_r < radius/earth_with_atm_radius) {
              return x_ye_max;
            }
            else if ( rel_r > radius/earth_with_atm_radius ){
              return 0.494;
            }else {
              return gsl_spline_eval(inter_ye,rel_r,inter_ye_accel);
            }
        }

EarthAtm::EarthAtm(std::string filepath):Body(7,"EarthAtm")
        {
            radius = 6371.0; // km
            atm_height = 100; // km
            earth_with_atm_radius = radius + atm_height;

            marray<double,2> earth_model = quickread(filepath);
            size_t arraysize = earth_model.extent(0);

            double earth_radius[arraysize];
            double earth_density[arraysize];
            double earth_ye[arraysize];

            for (unsigned int i=0; i < arraysize;i++){
                earth_radius[i] = earth_model[i][0];
                earth_density[i] = earth_model[i][1];
                earth_ye[i] = earth_model[i][2];
            }

            x_radius_min = earth_radius[0];
            x_radius_max = earth_radius[arraysize-1];
            x_rho_min = earth_density[0];
            x_rho_max = earth_density[arraysize-1];
            x_ye_min = earth_ye[0];
            x_ye_max = earth_ye[arraysize-1];

            inter_density = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_density_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_density,earth_radius,earth_density,arraysize);

            inter_ye = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_ye_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_ye,earth_radius,earth_ye,arraysize);
        }

EarthAtm::~EarthAtm(void)
        {
            gsl_spline_free(inter_density);
            gsl_interp_accel_free(inter_density_accel);
            gsl_spline_free(inter_ye);
            gsl_interp_accel_free(inter_ye_accel);
        }


}
