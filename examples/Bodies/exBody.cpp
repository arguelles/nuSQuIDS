#include "exBody.h"

#define SQR(x) ((x)*(x))

namespace nusquids{


static squids::Const param;

/*
----------------------------------------------------------------------
         EarthMod CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
  EarthMod::EarthMod():EarthMod(static_cast<std::string>(EARTH_MODEL_LOCATION),1,1,1)
  {
  }


  EarthMod::Track::Track(double phi_input):Body::Track(0,0)
  {
    radius_nu = 6371.0*param.km;
    //radius_nu = 6369.0*param.km;
    //atmheight = 100.0*param.km;
    atmheight = 22.*param.km;

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

#ifdef EarthMod_DEBUG
    cout << "== Init Track ==" << endl;
    cout << " phi = " << phi <<
      ", cos(phi) = " << cosphi <<
      ", L = " << radius_nu/param.km << endl;
    cout << "==" << endl;
#endif

    // TrackParams = {xini,xend,phi_input};
  }

  double EarthMod::density(const GenericTrack& track_input) const
  {
    const EarthMod::Track& track_earthatm = static_cast<const EarthMod::Track&>(track_input);
    double xkm = track_earthatm.GetX()/param.km;
    
    double r = sqrt(SQR(earth_with_atm_radius) + SQR(xkm) - (track_earthatm.L/param.km)*xkm);
    
#ifdef true
    //    cout << "************************************** r : " << r << " L : " << (track_earthatm->L/param.km)
    //	 << " x : " << xkm << " R : " << radius << endl;
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
  double EarthMod::ye(const GenericTrack& track_input) const
  {
    const EarthMod::Track& track_earthatm = static_cast<const EarthMod::Track&>(track_input);
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

  void EarthMod::Mod(double frho1, double frho2, double frho3){
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

  EarthMod::EarthMod(std::string filepath, double frho1, double frho2, double frho3):Body(10,"EarthMod")
  {
    // The Input file should have the radius specified from 0 to 1.
    // where 0 is the center of the EarthMod and 1 is the surface.
    atm_height = 20; // km
    radius = 6371.0; // [km]
    earth_with_atm_radius = radius + atm_height;

    earth_model = quickread(filepath);
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
  
  EarthMod::~EarthMod(){
    gsl_spline_free(inter_density);
    gsl_interp_accel_free(inter_density_accel);
    gsl_spline_free(inter_ye);
    gsl_interp_accel_free(inter_ye_accel);
  }
  

}
