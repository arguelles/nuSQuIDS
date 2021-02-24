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

#include <nuSQuIDS/body.h>

namespace nusquids {

// Exponential source
// The assumption will be that it will be a single flavor neutrino source
// with an exponential decay profile.
class EmittingVacuum: public Vacuum {
  private:
    const unsigned int flavor;
    const double decay_length;
  public :
    EmittingVacuum(unsigned int flavor, double decay_length):
      flavor(flavor),decay_length(decay_length){}
    void injected_neutrino_flux(marray<double,3>& flux, const GenericTrack& track, const nuSQUIDS& nusquids) override {
      double x_cur = track.GetX();
      for(unsigned int ei=0; ei < nusquids.GetNumE(); ei++){
        for(unsigned int rhoi = 0; rhoi < nusquids.GetNumRho(); rhoi++){
          for(unsigned int flv = 0; flv < nusquids.GetNumNeu(); flv++){
            flux[ei][rhoi][flv] = (flv == flavor) ? exp(-x_cur/decay_length) : 0.0;
          }
        }
      }
    }
  };
} // close nusquids namespace

#endif
