#include "exBody.h"

#define SQR(x) ((x)*(x))

namespace nusquids{


static squids::Const param;

/*
----------------------------------------------------------------------
         EarthMod CLASS DEFINITIONS
----------------------------------------------------------------------
*/


  void EarthMod::Mod(double frho1, double frho2, double frho3){
    marray<double,2> earth_model = quickread(static_cast<std::string>(EARTH_MODEL_LOCATION));
    size_t arraysize = earth_model.extent(0);
    double earth_radius[arraysize];
    double earth_density[arraysize];
    double earth_ye[arraysize];
    //std::ofstream file("model.txt", std::ofstream::app);
    //file << "# " << frho1 <<" "<< frho2 <<" " << frho3 << std::endl;
    for (unsigned int i=0; i < arraysize;i++){
      double d=1;
      if(i<39)
        d=frho1;
      else if (i>=39 && i<110)
        d=frho2;
      else if (i>=110 && i<201)
        d=frho3;
      
      earth_radius[i] = earth_model[i][0];
      earth_density[i] = d*earth_model[i][1];
      earth_ye[i] = earth_model[i][2];

      //  file << i <<" " << earth_radius[i] << " " << earth_density[i] << std::endl;  
    }
    //    file.close();

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

  EarthMod::EarthMod(std::string filepath, double frho1, double frho2, double frho3)
  {
    // The Input file should have the radius specified from 0 to 1.
    // where 0 is the center of the EarthMod and 1 is the surface.
    atm_height = 20; // km
    radius = 6371.0; // [km]
    earth_with_atm_radius = radius + atm_height;

    marray<double,2> earth_model = quickread(filepath);
    size_t arraysize = earth_model.extent(0);
    
    double earth_radius[arraysize];
    double earth_density[arraysize];
    double earth_ye[arraysize];
    
    for (unsigned int i=0; i < arraysize;i++){
      double d=1;
      if(i<39)
        d=frho1;
      else if (i>=39 && i<110)
        d=frho2;
      else if (i>=110 && i<201)
        d=frho3;

      earth_radius[i] = earth_model[i][0];
      earth_density[i] = d*earth_model[i][1];
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
  
}
