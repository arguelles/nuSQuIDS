#include <iostream>
#include <iomanip>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/tools.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using namespace nusquids;

int main(){
  // units
  const unsigned int numneu = 3;
  squids::Const units;
  double Emin = 1.0e2*units.GeV;
  double Emax = 1.0e6*units.GeV;
  unsigned int Esize = 100;
  nuSQUIDS nus(logspace(Emin,Emax,Esize),numneu,both,true);

  auto energies = nus.GetERange();
  nus.InitializeInteractions();
  auto is = nus.GetInteractionStructure();

  unsigned int nTargets = is->targets.size();
  for(unsigned int trg = 0; trg < nTargets; trg++){
  for(unsigned int rho = 0; rho < 2; rho++){
    for(unsigned int flv = 0; flv < numneu; flv++){
      for(unsigned int ei = 1; ei < nus.GetNumE(); ei++){
        // initialze muon state
        double total_CC_xs = is->sigma_CC[trg][rho][flv][ei];
        if(total_CC_xs < 0.)
          std::cout << "Error: negative cross sections." << std::endl;
        double total_NC_xs = is->sigma_NC[trg][rho][flv][ei];
        if(total_NC_xs < 0.)
          std::cout << "Error: negative cross sections." << std::endl;
        /*
        double total_dNdECC = 0;
        double begincap_dNdECC = energies[0]*is->dNdE_CC[rho][flv][ei][0];
        for(unsigned int ef = 0; ef < ei; ef++){
          total_dNdECC += (is->dNdE_CC[rho][flv][ei][ef+1]+is->dNdE_CC[rho][flv][ei][ef])*(energies[ef+1]-energies[ef])/2.;
        }
        if(std::abs(total_dNdECC + begincap_dNdECC - 1.0) > 1.0e-2)
          std::cout << "Error: CC cross sectiosn are not unitarity: "<< ei << " " << total_dNdECC + begincap_dNdECC << std::endl;
        */
        double total_dNdENC = 0;
        double begincap_dNdENC = energies[0]*is->dNdE_NC[trg][rho][flv][ei][0];
        double m = 0;
        for(unsigned int ef = 0; ef < ei; ef++){
          //total_dNdENC += is->dNdE_NC[rho][flv][ei][ef]*(energies[ef+1]-energies[ef]);
          //total_dNdENC += (is->dNdE_NC[rho][flv][ei][ef+1]+is->dNdE_NC[rho][flv][ei][ef])*(energies[ef+1]-energies[ef])/2.;
          if(ef + 1 == ei)
            total_dNdENC += (is->dNdE_NC[trg][rho][flv][ei][ef])*(energies[ef+1]-energies[ef]) + m*(energies[ef+1]-energies[ef])*(energies[ef+1]-energies[ef]);
          else
            total_dNdENC += (is->dNdE_NC[trg][rho][flv][ei][ef+1] + is->dNdE_NC[trg][rho][flv][ei][ef])*(energies[ef+1]-energies[ef])/2.;
          m = (is->dNdE_NC[trg][rho][flv][ei][ef+1]-is->dNdE_NC[trg][rho][flv][ei][ef])/(energies[ef+1]-energies[ef]);
          //std::cout << ef << " " << std::endl;
          //total_dNdENC += (is->dNdE_NC[rho][flv][ei][ef])*(energies[ef+1]-energies[ef]);
        }
        if(std::abs(total_dNdENC + begincap_dNdENC - 1.0) > 1.0e-2)
          std::cout << "Error: NC cross sectiosn are not unitarity: "<< rho << " " << flv << " " << ei << " " << total_dNdENC + begincap_dNdENC << std::endl;
      }
    }
  }
  }

  return 0;
}
