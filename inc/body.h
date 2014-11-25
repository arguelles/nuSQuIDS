#ifndef __BODY_H
#define __BODY_H

#include <string>
#include <math.h>
#include <cmath>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <SQuIDS/const.h>
#include "tools.h"
#include "global.h"
#include <assert.h>
#include <memory>
#include <exception>

namespace nusquids{

class Body{
    public:
        std::string name;
        int id;

        Body();

        class Track{
            public:
                double x;
                double xini;
                double xend;
                Track(double,double);
                Track();
                void SetX(double y){ x = y; };
                double GetX(){ return x; };
                double GetInitialX(){ return xini; };
                double GetFinalX(){ return xend; };
                std::vector<double> GetTrackParams(){ return TrackParams; };
                std::vector<double> TrackParams;
        };
        virtual double density(std::shared_ptr<Track>);
        virtual double ye(std::shared_ptr<Track>);
        const std::vector<double>& GetBodyParams(){ return BodyParams;};
        std::vector<double> BodyParams;
};

// type defining
typedef Body::Track GenericTrack;

class Vacuum: public Body {
    public:
        Vacuum();

        class Track: public Body::Track {
            public :
                Track(double,double);
                Track();
        };

        double density(std::shared_ptr<GenericTrack>);
        double ye(std::shared_ptr<GenericTrack>);
};

class ConstantDensity: public Body{
    public:
        double constant_density;
        double constant_ye;

        ConstantDensity(double,double);
        ConstantDensity();

        class Track: public Body::Track{
            public :
                Track(double,double);
                Track();
        };

        double density(std::shared_ptr<GenericTrack>);
        double ye(std::shared_ptr<GenericTrack>);
};

class VariableDensity: public Body{
    public:
        double* x_arr;
        double* density_arr;
        double* ye_arr;

        int arraysize;

        double x_max;
        double x_min;

        gsl_spline * inter_density;
        gsl_interp_accel * inter_density_accel;

        gsl_spline * inter_ye;
        gsl_interp_accel * inter_ye_accel;

        VariableDensity(std::vector<double>,std::vector<double>,std::vector<double>);
        VariableDensity();

        class Track: public Body::Track{
            public :
                Track(double,double);
                Track();
        };

        double density(std::shared_ptr<GenericTrack>);
        double ye(std::shared_ptr<GenericTrack>);
};

class Earth: public Body{
    public:
        double radius;
        double ye_mantle,ye_outercore,ye_innercore;

        gsl_spline * inter_density;
        gsl_interp_accel * inter_density_accel;
        gsl_spline * inter_ye;
        gsl_interp_accel * inter_ye_accel;

        double x_radius_min, x_radius_max;
        double x_rho_min, x_rho_max;
        double x_ye_min, x_ye_max;

        Earth();
        Earth(std::string);

        double rdensity(double);

        class Track: public Body::Track{
            public :
                double baseline;
                Track(double,double,double);
                Track();
        };

        double density(std::shared_ptr<Body::Track>);
        double ye(std::shared_ptr<GenericTrack>);

        ~Earth();
};

class Sun: public Body{
    public:
        Table sun_model;
        Table sun_model_nele;
        double* sun_radius;
        double* sun_density;
        double* sun_xh;
        double* sun_nele_radius;
        double* sun_nele;

        int arraysize,arraysize_2;

        double radius;

        gsl_spline * inter_density;
        gsl_interp_accel * inter_density_accel;

        gsl_spline * inter_rxh;
        gsl_interp_accel * inter_rxh_accel;

        gsl_spline * inter_nele;
        gsl_interp_accel * inter_nele_accel;

        Sun();

        double rdensity(double);
        double rxh(double);

        class Track: public Body::Track{
            public :
                Track(double,double);
                Track();
        };

        double density(std::shared_ptr<Body::Track>);
        double ye(std::shared_ptr<Body::Track>);
};

class SunASnu: public Body{
    public:
        Table sun_model;
        double* sun_radius;
        double* sun_density;
        double* sun_xh;

        int arraysize;

        double radius;

        gsl_spline * inter_density;
        gsl_interp_accel * inter_density_accel;

        gsl_spline * inter_rxh;
        gsl_interp_accel * inter_rxh_accel;

        SunASnu();

        double rdensity(double);
        double rxh(double);

        class Track: public Body::Track{
            public :
                double b_impact;
                double radius_nu;
                Track(double,double);
                Track();
        };

        double density(std::shared_ptr<Body::Track>);
        double ye(std::shared_ptr<Body::Track>);
};

class EarthAtm: public Body{
    public:
        double radius;
        double ye_mantle,ye_outercore,ye_innercore;

        gsl_spline * inter_density;
        gsl_interp_accel * inter_density_accel;
        gsl_spline * inter_ye;
        gsl_interp_accel * inter_ye_accel;

        double x_radius_min, x_radius_max;
        double x_rho_min, x_rho_max;
        double x_ye_min, x_ye_max;

        EarthAtm();
        EarthAtm(std::string);

        class Track: public Body::Track{
            public :
                double phi;
                double cosphi;
                double radius_nu;
                double atmheight;
                double L;

                Track(double);
                Track();
        };

        double density(std::shared_ptr<Body::Track>);
        double ye(std::shared_ptr<Body::Track>);
        ~EarthAtm();
};

// type defining
typedef Body::Track Track;

} // close namespace

#endif
