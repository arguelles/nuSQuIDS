#include <iostream>
#include "nuSQuIDS/nuSQuIDS.h"

using namespace nusquids;

int main(){
	squids::Const units;
	
	const unsigned int numneu=3;
	double Emin=1.e1*units.GeV;
	double Emax=1.e2*units.GeV;
	double czmin=-1;
	double czmax=0;
	nuSQUIDSAtm<> nus_atm(linspace(czmin,czmax,5),logspace(Emin,Emax,10),numneu,both,false);
	
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
	nus_atm.Set_IncludeOscillations(true);
	
	nus_atm.EvolveState();
	
	//simple sanity check: all fluxes should remain in [0,1] for this simple case
	for(double cz : linspace(czmin,czmax,20)){
		for(double en : logspace(Emin,Emax,100)){
			for(int fl=0; fl<numneu; fl++){
				double flux=nus_atm.EvalFlavor(fl,cz,en);
				if(flux<0 || flux>1)
					std::cout << "Bad flux: " << flux << " for cos(angle)=" <<
					cz << ", energy=" << en/units.GeV << " GeV, flavor " << fl
					<< std::endl;
			}
		}
	}
}