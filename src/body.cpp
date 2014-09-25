#include "body.h"

// Macros
#define SQR(x)      ((x)*(x))                        // x^2
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  // |x|^2
#define POW10(x)    (exp(M_LN10*(x)))                // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( (b) > 0.0 ? (fabs(a)) : (-fabs(a)) )
#define KRONECKER(i,j)  ( (i)==(j) ? 1 : 0 )

Const param;

/*
-----------------------------------------------------------------------
         BODY CLASS DEFINITIONS
-----------------------------------------------------------------------
*/


Body::Body()
        {
            id = 0;
            name = "GenericBody";
        }

// track constructor
Body::Track::Track(double xini_input, double xend_input): TrackParams({xini_input,xend_input})
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
        }
Body::Track::Track(): TrackParams({})
        {
            //pass
        }

double Body::density(std::shared_ptr<Track> track_input){
            return 0.0;
        }

double Body::ye(std::shared_ptr<Track> track_input){
            return 1.0;
        }

/*
----------------------------------------------------------------------
         VACUUM CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// vacuum constructor
Vacuum::Vacuum()
        {
            id = 1;
            name = "Vacuum";
        }

// track constructor
Vacuum::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
            TrackParams = {xini,xend};
        }
Vacuum::Track::Track()
        {
            //pass
        }

double Vacuum::density(std::shared_ptr<GenericTrack> track_input){
            return 0.0;
        }

double Vacuum::ye(std::shared_ptr<GenericTrack> track_input){
            return 1.0;
        }

/*
----------------------------------------------------------------------
         ConstantDensity CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
ConstantDensity::ConstantDensity(double density_input,double ye_input)
        {
            name = "ConstantDensity";
            id = 2;
            constant_density = density_input;
            constant_ye = ye_input;
            BodyParams = {constant_density, constant_ye};
        }

ConstantDensity::ConstantDensity(void)
        {
        }

// track constructor
ConstantDensity::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
            TrackParams = {xini,xend};
        }
ConstantDensity::Track::Track()
        {
            //pass
        }

double ConstantDensity::density(std::shared_ptr<GenericTrack> track_input)
        {
            return constant_density;
        }

double ConstantDensity::ye(std::shared_ptr<GenericTrack> track_input)
        {
            return constant_ye;
        }

/*
----------------------------------------------------------------------
         VariableDensity CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
VariableDensity::VariableDensity(std::vector<double> x_input,std::vector<double> density_input,std::vector<double> ye_input)
        {
            name = "VariableDensity";
            id = 3;

            assert("nuSQUIDS::Error::VariableDensityConstructor: Invalid array sizes." && x_input.size() == density_input.size() && x_input.size() == ye_input.size());
            arraysize = x_input.size();

            x_min = x_input.front();
            x_max = x_input.back();

            x_arr = new double[arraysize];
            density_arr = new double[arraysize];
            ye_arr = new double [arraysize];

            for(int i = 0; i < arraysize; i++){
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

VariableDensity::VariableDensity(){};
// track constructor
VariableDensity::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
            TrackParams = {xini,xend};
        }
VariableDensity::Track::Track()
        {
            //pass
        }

double VariableDensity::density(std::shared_ptr<GenericTrack> track_input)
        {
          double x = track_input->x;
          if (x < x_min or x > x_max ){
              return 0;
          } else {
              return gsl_spline_eval(inter_density,x,inter_density_accel);
          }
        }
double VariableDensity::ye(std::shared_ptr<GenericTrack> track_input)
        {
          double x = track_input->x;
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
Earth::Earth()
        {
            name = "Earth";
            id = 4;
            radius = 6371.0; // [km]
            ye_mantle = 0.4957;
            ye_outercore = 0.4656;
            ye_innercore = 0.4656;
        }

/*
Earth::Earth(string earth_table_path)
        {


        }
*/
// track constructor
Earth::Track::Track(double xini_input, double xend_input,double baseline_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
            baseline = baseline_input;
            TrackParams = {xini,xend,baseline};
        }
Earth::Track::Track()
        {
            //pass
        }

double Earth::rdensity(double x){
        // Calcula la densidad de la Tierra segun el PREM
        // R. Gandhi et al. Astroparticle Phys. 5, 81-110 (1996)
        // Arxiv : 9512364 pag. 23
        // x is adimentional radius : x = 0 : center, x = 1 : Earthradius
        double dne = 0.0;
        double r = radius*x;
        if (r <= 1221.50)
            dne = 13.08850 - 8.83810*SQR(x);
        else if (r>=1221.50 and r<3480) 
            dne=12.58150 - 1.26380*x - 3.64260*SQR(x) - 5.5280*SQR(x)*x;
        else if (r >=3480.0 and r < 5701.0)
            dne=7.95650 - 6.47610*x + 5.52830*SQR(x) - 3.08070*SQR(x)*x;
        else if (r >= 5701.0 and r<5771.0)
            dne=5.31970 - 1.48360*x;
        else if (r>=5771.0 and r<5971.0)
            dne=11.24940 - 8.02980*x;
        else if (r>=5971.0 and r<6151.0)
            dne=7.10890 - 3.80450*x;
        else if (r>=6151.0 and r<6346.60)
            dne=2.6910 + 0.69240*x;
        else if (r >= 6346.60 and r < 6356.0)
            dne = 2.9;
        else if (r >= 6356.0 and r < 6368)
            dne = 2.6;
        else if (r<= radius)
            dne = 1.020;
        else if (r>= radius)
            dne=0.0;
        return dne;
        }

double Earth::density(std::shared_ptr<Body::Track> track_input)
        {
            std::shared_ptr<Earth::Track> track_earth = std::static_pointer_cast< Earth::Track >(track_input);
            double xkm = track_earth->x/param.km;
            double r = sqrt(SQR(radius)+SQR(xkm)-(track_earth->baseline/param.km)*xkm);
            if (use_file){
              //if ( r/radius < x_radius_min or r/radius > x_radius_max)
              //  return 0.0;
              return gsl_spline_eval(inter_density,r/radius,inter_density_accel);
            }
            else{
              return rdensity(r/radius);
            }
        }
double Earth::ye(std::shared_ptr<Body::Track> track_input)
        {
            std::shared_ptr<EarthAtm::Track> track_earthatm = std::static_pointer_cast< EarthAtm::Track >(track_input);
            if(track_earthatm->cosphi <= 0.0){
                double xkm = track_earthatm->x/param.km;
                double r = sqrt(SQR(radius) + SQR(xkm) - (track_earthatm->L/param.km)*xkm);

                if (use_file){
                  //if ( r/radius < x_radius_min or r/radius > x_radius_max)
                  //  return 0.0;
                  return gsl_spline_eval(inter_ye,r/radius,inter_ye_accel);
                }

                if ( r < 1121.5)
                  return ye_innercore;
                else if ( r < 3480. )
                  return ye_outercore;
                else
                  return ye_mantle;
            } else {
              return 0.494;
            }
        }

Earth::Earth(string filepath)
        {
          // The Input file should have the radius specified from 0 to 1.
          // where 0 is the center of the Earth and 1 is the surface.
            use_file = true;

            name = "Earth";
            id = 4;
            radius = 6371.0; // [km]

            Table earth_model = quickread(filepath);
            size_t arraysize = earth_model.size();

            double earth_radius[arraysize];
            double earth_density[arraysize];
            double earth_ye[arraysize];

            for (int i=0; i < arraysize;i++){
                earth_radius[i] = earth_model[i][0];
                earth_density[i] = earth_model[i][1];
                earth_ye[i] = earth_model[i][2];
            }

            double x_radius_min = earth_radius[0];
            double x_radius_max = earth_radius[arraysize];

            inter_density = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_density_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_density,earth_radius,earth_density,arraysize);

            inter_ye = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_ye_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_ye,earth_radius,earth_ye,arraysize);
        }

Earth::~Earth(){
  if(use_file){
    free(inter_density);
    free(inter_density_accel);
    free(inter_ye);
    free(inter_ye_accel);
  }
}

/*
----------------------------------------------------------------------
         SUN CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
Sun::Sun()
        {
            name = "Sun";
            id = 5;
            //radius = 694439.0;
            radius = 695980.0;


            // import sun model
            sun_model = quickread(SUN_MODEL_LOCATION);
            arraysize = sun_model.size();

            sun_radius = new double[arraysize];
            sun_density = new double[arraysize];
            sun_xh = new double[arraysize];

            for (int i=0; i < arraysize;i++){
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
Sun::Track::Track(double xini_input, double xend_input)
        {
            x = xini_input;
            xini = xini_input;
            xend = xend_input;
            TrackParams = {xini,xend};
        }
Sun::Track::Track()
        {
            //pass
        }
double Sun::rdensity(double x){
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_density[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_density,x,inter_density_accel);
            }
        }

double Sun::rxh(double x){
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_xh[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_rxh,x,inter_rxh_accel);
            }
        }

double Sun::density(std::shared_ptr<GenericTrack> track_input)
        {
            double r = track_input->x/(radius*param.km);
            return rdensity(r);
        }
double Sun::ye(std::shared_ptr<GenericTrack> track_input)
        {
            double r = track_input->x/(radius*param.km);
            return 0.5*(1.0+rxh(r));
        }

/*
----------------------------------------------------------------------
         SUN ASNU CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
SunASnu::SunASnu()
        {
            name = "SunASnu";
            id = 6;
            radius = 694439.0;

            sun_model = quickread(SUN_MODEL_LOCATION);
            arraysize = sun_model.size();

            sun_radius = new double[arraysize];
            sun_density = new double[arraysize];
            sun_xh = new double[arraysize];

            for (int i=0; i < arraysize;i++){
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
SunASnu::Track::Track(double xini_input, double b_impact_input)
        {
            radius_nu = 694439.0*param.km;
            x = xini_input;
            xini = xini_input;
            b_impact = b_impact_input;
            xend = 2.0*sqrt(SQR(radius_nu)+SQR(b_impact));
            TrackParams = {xini,xend, b_impact_input};
        }
SunASnu::Track::Track()
        {
            //pass
        }

double SunASnu::rdensity(double x){
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_density[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_density,x,inter_density_accel);
            }
        }

double SunASnu::rxh(double x){
        // x is adimentional radius : x = 0 : center, x = 1 : radius
            if (x < sun_radius[0]){
                return sun_xh[0];
            } else if ( x > sun_radius[arraysize-1]){
                return 0;
            } else {
                return gsl_spline_eval(inter_rxh,x,inter_rxh_accel);
            }
        }

double SunASnu::density(std::shared_ptr<Body::Track> track_input)
        {
            std::shared_ptr<SunASnu::Track> track_sunasnu = std::static_pointer_cast< SunASnu::Track >(track_input);
            double xkm = track_sunasnu->x/param.km;
            double bkm = track_sunasnu->b_impact/param.km;

            double r = sqrt(SQR(radius)+SQR(xkm)-2.0*xkm*sqrt(SQR(radius)-SQR(bkm)))/radius;

            return rdensity(r);
        }

double SunASnu::ye(std::shared_ptr<Body::Track> track_input)
        {
            std::shared_ptr<SunASnu::Track> track_sunasnu = std::static_pointer_cast< SunASnu::Track >(track_input);
            double xkm = track_sunasnu->x/param.km;
            double bkm = track_sunasnu->b_impact/param.km;
            double r = sqrt(SQR(radius)+SQR(xkm)-2.0*xkm*sqrt(SQR(radius)-SQR(bkm)))/radius;
            return 0.5*(1.0+rxh(r));
        }

/*
----------------------------------------------------------------------
         EARTHATM CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
EarthAtm::EarthAtm(void)
        {
            name = "EarthAtm";
            id = 7;
            radius = 6371.0;

            ye_mantle = 0.4957;
            ye_outercore = 0.4656;
            ye_innercore = 0.4656;
        }
// track constructor
EarthAtm::Track::Track(double phi_input)
        {
            radius_nu = 6371.0*param.km;
            atmheight = 100.0*param.km;

            phi = phi_input;
            cosphi = cos(phi);


            if(cosphi<=0.0){
                L = 2.0*radius_nu*abs(cosphi);
            } else {
                L = atmheight/abs(cosphi);
            }


            double R = radius_nu;
            double r = atmheight;
            double mm = tan(phi);

            /*
            if(cosphi<=0.0){
                L = sqrt(((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*r*R + 2.0*R*
                    (R + sqrt((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*
                    r*R + R*R)))/(1.0 + mm*mm));
            } else {
                L = sqrt(((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*r*R + 2.0*R*
                    (R - sqrt((1.0 + mm*mm)*r*r + 2.0*(1.0 + mm*mm)*
                    r*R + R*R)))/(1.0 + mm*mm));
            }
            */

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

            TrackParams = {xini,xend, phi_input};
        }
EarthAtm::Track::Track()
        {
            //pass
        }
double EarthAtm::rdensity(double x){
        // Calcula la densidad de la Tierra segun el PREM
        // R. Gandhi et al. Astroparticle Phys. 5, 81-110 (1996)
        // Arxiv : 9512364 pag. 23
        // x is adimentional radius : x = 0 : center, x = 1 : Earthradius
        double dne = 0.0;
        double r = radius*x;
        if (r <= 1221.50)
            dne = 13.08850 - 8.83810*SQR(x);
        else if (r>=1221.50 and r<3480) 
            dne=12.58150 - 1.26380*x - 3.64260*SQR(x) - 5.5280*SQR(x)*x;
        else if (r >=3480.0 and r < 5701.0)
            dne=7.95650 - 6.47610*x + 5.52830*SQR(x) - 3.08070*SQR(x)*x;
        else if (r >= 5701.0 and r<5771.0)
            dne=5.31970 - 1.48360*x;
        else if (r>=5771.0 and r<5971.0)
            dne=11.24940 - 8.02980*x;
        else if (r>=5971.0 and r<6151.0)
            dne=7.10890 - 3.80450*x;
        else if (r>=6151.0 and r<6346.60)
            dne=2.6910 + 0.69240*x;
        else if (r >= 6346.60 and r < 6356.0)
            dne = 2.9;
        else if (r >= 6356.0 and r < 6368.0)
            dne = 2.6;
        else if (r<= radius)
            dne = 1.020;
        else if (r>= radius)
            dne=0.0;
        return dne;
        }

double EarthAtm::density(std::shared_ptr<Body::Track> track_input)
        {
            std::shared_ptr<EarthAtm::Track> track_earthatm = std::static_pointer_cast< EarthAtm::Track >(track_input);
            double xkm = track_earthatm->x/param.km;

            if(track_earthatm->cosphi <= 0.0){
                double r = sqrt(SQR(radius) + SQR(xkm) - (track_earthatm->L/param.km)*xkm);

                #ifdef EarthAtm_DEBUG
                cout << "r : " << r << " L : " << (track_earthatm->L/param.km)
                     << " x : " << xkm << " R : " << radius << endl;
                #endif

                if (use_file){
                  // fix latter
                  //if ( r/radius < x_radius_min or r/radius > x_radius_max)
                  //  return 0.0;
                  return gsl_spline_eval(inter_density,r/radius,inter_density_accel);
                }

                return rdensity(r/radius);
            } else {
                double h = xkm*track_earthatm->cosphi;
                double h0 = 25.0;
                return 1.05*exp(h/h0);
            }
        }
double EarthAtm::ye(std::shared_ptr<Body::Track> track_input)
        {
            std::shared_ptr<EarthAtm::Track> track_earthatm = std::static_pointer_cast< EarthAtm::Track >(track_input);
            if(track_earthatm->cosphi <= 0.0){
                double xkm = track_earthatm->x/param.km;
                double r = sqrt(SQR(radius) + SQR(xkm) - (track_earthatm->L/param.km)*xkm);

                if (use_file){
                  //if ( r/radius < x_radius_min )
                  //  or r/radius > x_radius_max)
                  //  return 0.0;
                  return gsl_spline_eval(inter_ye,r/radius,inter_ye_accel);
                }

                if ( r < 1121.5)
                  return ye_innercore;
                else if ( r < 3480. )
                  return ye_outercore;
                else
                  return ye_mantle;
            } else {
              return 0.494;
            }
        }

EarthAtm::EarthAtm(string filepath)
        {
            use_file = true;

            name = "EarthAtm";
            id = 7;
            radius = 6371.0;

            Table earth_model = quickread(filepath);
            size_t arraysize = earth_model.size();

            double earth_radius[arraysize];
            double earth_density[arraysize];
            double earth_ye[arraysize];

            for (int i=0; i < arraysize;i++){
                earth_radius[i] = earth_model[i][0];
                earth_density[i] = earth_model[i][1];
                earth_ye[i] = earth_model[i][2];
            }

            inter_density = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_density_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_density,earth_radius,earth_density,arraysize);

            inter_ye = gsl_spline_alloc(gsl_interp_cspline,arraysize);
            inter_ye_accel = gsl_interp_accel_alloc ();
            gsl_spline_init (inter_ye,earth_radius,earth_ye,arraysize);
        }

EarthAtm::~EarthAtm(void)
        {
            if(use_file){
              free(inter_density);
              free(inter_density_accel);
              free(inter_ye);
              free(inter_ye_accel);
            }
        }
