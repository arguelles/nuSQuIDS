#ifndef _SMinit_
#define _SMinit_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include <memory>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "marray.h"
#include "tools.h"

class littlemermaid {
  public:
    /// \brief component settings
    enum FluxType {pp = 0, pep = 1, hep = 2, be7 = 3, b8 = 4, n13 = 5, o15 = 6, f17 = 7, Electron = 8, DM = 9};
    /// \brief line features
    std::map<FluxType,bool> isline {{pp,false},{pep,true},{hep,false},{be7,true},{n13,false},{o15,false},{f17,false}};
	private:
    const std::string datapath = "SM/";
    const double Emin = 3.4640e3;
    const double Emax = 1.8784e07;
    const double Rsun = 69662651411.1;
    const double MeV = 1.0e6;
    const unsigned int num_components = 10;
    std::map<FluxType,std::string> spectrum_filename {{pp,"pp.dat"},{pep,""},{hep,"hep.dat"},{be7,""},{b8,"b8.dat"},{n13,"n13.dat"},{o15,"o15.dat"},{f17,"f17.dat"},{Electron,""},{DM,""}};
  public:
    std::map<FluxType,std::vector<double>> spectrum_limits {{pp,{0.00504e6,0.42341e6}},
                                                            {pep,{1.44e6}},
                                                            {hep,{1.8784e4,1.8784e7}},
                                                            {be7,{0.8618e6}},
                                                            {b8,{0.02e6,16.56e6}},
                                                            {n13,{5.9950e3,1.1990e6}},
                                                            {o15,{3.4640e3,1.7320e6}},
                                                            {f17,{3.4800e3,1.7400e6}},
                                                            {Electron,{}},
                                                            {DM,{}}};
  private:
    ///\brief Splines
    std::vector<std::shared_ptr<gsl_interp_accel>> Racc;
    std::vector<std::shared_ptr<gsl_spline>> Rspline;

    std::vector<std::shared_ptr<gsl_interp_accel>> Eacc;
    std::vector<std::shared_ptr<gsl_spline>> Espline;
    std::string SM;
    ///\brief Spline initializer helper function
		void splineinit();
	public:
		littlemermaid(std::string solarmodel):SM(solarmodel){
      Racc.resize(num_components);Rspline.resize(num_components);
      Eacc.resize(num_components);Espline.resize(num_components);
      for(auto& racc : Racc){
        racc.reset(gsl_interp_accel_alloc(),[](gsl_interp_accel* t){ gsl_interp_accel_free(t);});
      }
      for(auto& eacc : Eacc){
        eacc.reset(gsl_interp_accel_alloc(),[](gsl_interp_accel* t){ gsl_interp_accel_free(t);});
      }

      // initialize splines
			splineinit();
		}
		double nuFlux(double r, double ee, FluxType) const;
		double eDensity(double r) const;
		double DMDensity(double r) const;
    unsigned int NumComp() const;

    //virtual ~littlemermaid();
};

#endif
