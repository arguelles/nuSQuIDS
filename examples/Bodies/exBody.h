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
#include <nuSQuIDS/resources.h>

/* This is a simple example where we modified the earth object, constructing a new object
 * called EarthMod, this very simple modification consist in including 3 parameters that weight 
 * density in every layer of the standard earth model.
 */

namespace nusquids {

class EarthMod: public EarthAtm{
public:
  //Constructor rescaling the densities in every layer
  EarthMod(std::string earthmodel, double frho1, double frho2 , double frho3);
  //Constructor rescaling the densities in every layer with default PREM.
  EarthMod(double frho1, double frho2 , double frho3):EarthMod(getResourcePath()+"/astro/EARTH_MODEL_PREM.dat",frho1,frho2,frho3){};
  //This function sets the values that weight the different layers in the of the PREM model
  void  Mod(double frho1, double frho2, double frho3);
};

} // close nusquids namespace

#endif
