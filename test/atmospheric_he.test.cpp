#include <iostream>
#include "nuSQuIDS/nuSQuIDS.h"

using namespace nusquids;

int main(){
	squids::Const units;
	
	const unsigned int numneu=3;
	double Emin=1.e5*units.GeV;
	double Emax=1.e7*units.GeV;
	double czmin=-1;
	double czmax=0;
	nuSQUIDSAtm<> nus_atm(linspace(czmin,czmax,5),logspace(Emin,Emax,10),numneu,both,true);
	
	nus_atm.Set_MixingAngle(0,1,0.563942);
	nus_atm.Set_MixingAngle(0,2,0.154085);
	nus_atm.Set_MixingAngle(1,2,0.785398);
	nus_atm.Set_SquareMassDifference(1,7.65e-05);
	nus_atm.Set_SquareMassDifference(2,0.00247);
	nus_atm.Set_CPPhase(0,2,0);
	
	marray<double,4> inistate{nus_atm.GetNumCos(),nus_atm.GetNumE(),2,numneu};
	std::fill(inistate.begin(),inistate.end(),1); //unit flux at all energies and angles
	nus_atm.Set_initial_state(inistate,flavor);
	nus_atm.Set_ProgressBar(false);
	nus_atm.Set_IncludeOscillations(false);
	nus_atm.Set_TauRegeneration(true);
	nus_atm.Set_GlashowResonance(true);
	
	nus_atm.Set_rel_error(1.0e-10);
	nus_atm.Set_abs_error(1.0e-10);
	
	nus_atm.EvolveState();
	
	//very simple sanity check: all fluxes should be real and effectively non-negative
	for(double cz : linspace(czmin,czmax,20)){
		for(double en : logspace(Emin,Emax,100)){
			for(int fl=0; fl<numneu; fl++){
				double flux=nus_atm.EvalFlavor(fl,cz,en);
				if(std::isnan(flux) || flux<-1e-4)
					std::cout << "Bad flux: " << flux << " for cos(angle)=" <<
					cz << ", energy=" << en/units.GeV << " GeV, flavor " << fl
					<< std::endl;
			}
		}
	}
}