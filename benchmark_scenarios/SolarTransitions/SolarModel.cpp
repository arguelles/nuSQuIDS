#include <SMinit.h>

using namespace std;
using namespace nusquids;

void littlemermaid::splineinit(){
	string rad = datapath + SM + "/radial.dat";

	marray<double,2> fluxR = quickread(rad);

	unsigned int arraysize = fluxR.extent(0);
	double XfluxRarr[arraysize];
  	for(unsigned int j = 1; j < arraysize; j++){
    	XfluxRarr[j] = fluxR[j][0];
  	}

	for(unsigned int i = 0; i <= 9; i++){
    	Rspline[i].reset(gsl_spline_alloc (gsl_interp_linear, arraysize),[](gsl_spline* t){gsl_spline_free(t);});
    	double YfluxRarr[arraysize];
    	for(unsigned int j = 0; j < arraysize; j++){
	        YfluxRarr[j] = fluxR[j][i+1];
	    }
    	gsl_spline_init (Rspline[i].get(), XfluxRarr, YfluxRarr, arraysize);
	}

	for(unsigned int i = 0; i < num_components; i++){
		if ( spectrum_filename[FluxType(i)] == ""){
		Espline[i] = nullptr;
		continue;
	}
	    marray<double,2> spectrum_table = quickread(datapath + "/NuSpec/" + spectrum_filename[FluxType(i)]);
	    unsigned int arraysize_ = spectrum_table.extent(0);
	    Espline[i].reset(gsl_spline_alloc (gsl_interp_linear, arraysize_),[](gsl_spline* t){gsl_spline_free(t);});
	    double ESpectrum[arraysize_]; double FSpectrum[arraysize_];
	    for(unsigned int j = 0; j < arraysize_ ; j++){
		      ESpectrum[j] = spectrum_table[j][0]*MeV;
		      FSpectrum[j] = spectrum_table[j][1];
	    }
	    gsl_spline_init(Espline[i].get(),ESpectrum,FSpectrum,arraysize_);
	}
}

double littlemermaid::nuFlux(double R, double E, FluxType type) const{
	if (R > 1.0 || R < 0.0)
    	throw std::runtime_error("Invalid R, must be between Rmin and Rmax");
	if (E < Emin || E > Emax){
		return 0.0;
  	}
  	if ( spectrum_limits.find(type)->second.size() == 1 ){
    	return gsl_spline_eval(Rspline[type].get(),R*Rsun,Racc[type].get());
 	} else if ( spectrum_limits.find(type)->second.size() == 2){
    if ( E > spectrum_limits.find(type)->second.front() and E < spectrum_limits.find(type)->second.back() )
    	return gsl_spline_eval(Rspline[type].get(),R*Rsun,Racc[type].get())*gsl_spline_eval(Espline[type].get(),E,Eacc[type].get());
    else
    	return 0.;
    } else{
    	throw std::runtime_error("Incorrect type, pp-0, pep-1, hep-2, be7-3, b8-4, N13-5, O15-6, F17-7");
    }
}

double littlemermaid::eDensity(double R) const {
	if (R > 1.0 || R < 0.0){
    throw std::runtime_error("Invalid R, must be between 0 and 1.");
  }
	return gsl_spline_eval(Rspline[Electron].get(),R*Rsun,Racc[Electron].get());
}

double littlemermaid::DMDensity(double R) const {
	if (R > 1.0 || R < 0.0){
		throw std::runtime_error("Invalid R, must be between 0 and 1.");
  }
	return gsl_spline_eval(Rspline[DM].get(),R*Rsun,Racc[DM].get());
}

unsigned int littlemermaid::NumComp() const {
  return 8;
}
