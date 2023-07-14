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
#include "NSI.h"
#include <nuSQuIDS/body.h>
#include <SQuIDS/const.h>
#include <nuSQuIDS/nuSQuIDS.h>

#define SQR(x)      ((x)*(x)) 

class EmittingEarthAtm: public EarthAtm {
    private:

    public:    
          EmittingEarthAtm( std::optional<std::string> atm_emission_filename = std::nullopt ):
          atm_emission_filename(atm_emission_filename) {}
      void injected_neutrino_flux(marray<double, 3>& flux, const GenericTrack& track, const nuSQUIDS& nusquids) override {


    std::vector<double> cz_data;
    std::vector<double> E_data;
    std::vector<double> h_data;

    
    std::ifstream inputFile2("czenfiles/" + czendata);
    double cz;
    double h;
    double E;

// add other flavors or parameters which your emission depends on here as they appear in your data
    while (inputFile2 >> cz >> h >> E >> nuEflux >> nuMuflux) {
      cz_data.push_back(cz);
      h_data.push_back(h);
      E_data.push_back(E);
    }

    inputFile2.close();

    int Unique_Sort(std::vector<double>& list) {
      std::sort(list.begin(), list.end());
      auto last = std::unique(lsit.begin(), list.end());
      list.erase(last, list.end());
    }
    // Sort the vectors
    Unique_Sort(cz_data.begin());
    Unique_Sort(h_data.begin());
    Unique_Sort(E_data.begin());
    double interpolation_min_height = h_data.front();
    double interpolation_max_height = h_data.front();
    double interpolation_min_Energy = E_data.back();
    double interpolation_max_Energy = E_data.back();

    marray<double, 3> nuEmatrices{cz_data.size(), h_data.size(), E_data.size()};
    marray<double, 3> nuMumatrices{cz_data.size(), h_data.size(), E_data.size()};

    std::ifstream inputFile2("czenfiles/" + czendata);
    double cz_;
    double h_;
    double E_;
    double nuEflux;
    double nuMuflux;

    while (inputFile2 >> cz_ >> h_ >> E_ >> nuEflux >> nuMuflux) {
        auto it1 = std::find(cz_data.begin(), cz_data.end(), cz_);
        auto it2 = std::find(h_data.begin(), h_data.end(), h_);
        auto it3 = std::find(E_data.begin(), E_data.end(), E_);
        if (it1 != cz_data.end() && it2 != h_data.end() && it3 != E_data.end()) {
            size_t czIndex = std::distance(cz_data.begin(), it1); 
            size_t hIndex = std::distance(h_data.begin(), it2);
            size_t EIndex = std::distance(E_data.begin(), it3);
            nuEmatrices[czIndex][hIndex][EIndex] = nuEflux;
            nuMumatrices[czIndex][hIndex][EIndex] = nuMuflux;
        }
    }
    inputFile2.close();
        
    marray<double> h_array(h_data.size());
    std::copy(h_data.begin(), h_data.end(), h_array.begin());
    marray<double> E_array(E_data.size());
    std::copy(E_data.begin(), E_data.end(), E_array.begin());
        
    marray<BicubicInterpolator, 2> Interpolators{cz_data.size, 2};

    for (int i = 0; i < cz_data.size, i++) {
      BiCubicInterpolator nuEinterpolator(nuEmatrices[i], E_array, h_array);
      BiCubicInterpolator nuMuinterpolator(nuMumatrices[i], E_array, h_array);
      Interpolators[i][0].push_back(nuEinterpolator);
      Interpolators[i][1].push_back(nuMuinterpolator);
    }
 
    // height calculation
    double xkm = track.GetX() / param.km;
    double sinsqphi = 1 - SQR(czNode);
    double dL = sqrt(SQR(earth_with_atm_radius) - SQR(radius) * sinsqphi) + radius * czNode;
    double L = sqrt(SQR(radius + atm_height) - SQR(radius) * sinsqphi) - radius * czNode;
    double r2 = SQR(earth_with_atm_radius) + SQR(xkm) - (L / param.km + dL) * xkm;
    double r = (r2 > 0 ? sqrt(r2) : 0);
    double rel_r = r / earth_with_atm_radius;
    double Height = earth_with_atm_radius * (rel_r - radius / earth_with_atm_radius); // here atm_height is a member of EarthAtm in km

    //setting the flux to the correct function, 0 outside range of interpolation
      if (Height >= interpolation_min_height && Height <= interpolation_max_height) {    
        for (unsigned int ei=0; ei < nusquids.GetNumE(); ei++) {
          double Energy = ESpace[ei];
          if (Energy >= interpolation_min_Energy && Height <= interpolation_max_Energy){
            for(unsigned int rhoi = 0; rhoi < nusquids.GetNumRho(); rhoi++){
              flux[ei][rhoi][4] = nuEinterpolator(Energy, Height*100000); //(height*100000 passed in cm)
              flux[ei][rhoi][5] = nuMuinterpolator(Energy, Height*100000); //(height*100000 passed in cm)
            }
          }
          else {
            for (unsigned int ei=0; ei < nusquids.GetNumE(); ei++) {
              for(unsigned int rhoi = 0; rhoi < nusquids.GetNumRho(); rhoi++){
                flux[ei][rhoi][4] = 0;
                flux[ei][rhoi][5] = 0;
              }
            }
          }
        }
      }
      else {
        for (unsigned int ei=0; ei < nusquids.GetNumE(); ei++) {
          for(unsigned int rhoi = 0; rhoi < nusquids.GetNumRho(); rhoi++){
            flux[ei][rhoi][4] = 0;
            flux[ei][rhoi][5] = 0;
          }
        }
      }

        
    }
  };


#endif
