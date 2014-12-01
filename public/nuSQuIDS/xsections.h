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

class NeutrinoCrossSections{
    private :
      double Emin;
      double Emax;
      int div;

      Table dsde_CC_tbl;
      Table dsde_NC_tbl;
      Table sigma_CC_tbl;
      Table sigma_NC_tbl;

      Table dsde_CC_data,dsde_NC_data,sigma_CC_data,sigma_NC_data;

      double LinInter(double,double,double,double,double);
      double SigmaInter(double, gsl_spline *, gsl_interp_accel *);
      double DSDEInter(double,double,int,std::vector<double>,std::string);
    public :
      NeutrinoCrossSections();
      NeutrinoCrossSections(double,double,int);
      void Init(double,double,int);

      double dsde_CC(int,int,int,int);
      double dsde_NC(int,int,int,int);
      double sigma_CC(int,int,int);
      double sigma_NC(int,int,int);
};

} // close namespace

#endif
