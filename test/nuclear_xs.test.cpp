#include <algorithm>
#include <cassert>
#include <iostream>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/resources.h>

int main(){
	using namespace nusquids;
	squids::Const units;
	
	auto csms=loadDefaultCrossSections();
	auto wcg24Nucleon=loadWCG24CrossSections();
	auto wcg24Nuclear=std::make_shared<CrossSectionLibrary>(loadWCG24NuclearCrossSections());
	
	//targets required to be able to approximate any normal matter
	auto generalTargets={electron,proton,neutron};
	{
		auto csmsTargets=csms.targets();
		for(auto target : generalTargets){
			if(!csmsTargets.count(target))
				std::cout << "CSMS cross section library is missing basic target " << target << '\n';
		}
		auto wcg24Targets=wcg24Nucleon.targets();
		for(auto target : generalTargets){
			if(!wcg24Targets.count(target))
				std::cout << "WCG24 cross section library is missing basic target " << target << '\n';
		}
		//wcg24Nuclear is not usable for arbitrary materials so we do not check it
	}
	
	{ //wcg24Nuclear should have some targets detected as nuclear (but not all, due to electrons)
		auto wcg24NuclearTargets=wcg24Nuclear->targets();
		assert(std::any_of(wcg24NuclearTargets.begin(), wcg24NuclearTargets.end(), isNuclearPDGCode));
	}
	
	{ //construct an EarthAtm with default (nucleon) cross sections, and ensure it does not use nuclear targets
		nuSQUIDS nus(marray<double, 1>({3},{100*units.GeV, 1e3*units.GeV, 1e4*units.GeV}), 3, neutrino, true);
		//leave cross sections as default (CSMS)
		std::string prem_path=getResourcePath()+"/astro/EARTH_MODEL_PREM.dat";
        auto earth=std::make_shared<EarthAtm>(prem_path);
        nus.Set_Body(earth);
        auto track=earth->MakeTrackWithCosine(-1);
        nus.Set_Track(std::make_shared<EarthAtm::Track>(track));
        nus.InitializeInteractions();
        auto intData=nus.GetInteractionStructure();
        assert(intData);
        assert(std::none_of(intData->targets.begin(), intData->targets.end(), isNuclearPDGCode));
	}
	
	{ //construct an EarthAtm with nuclear cross sections, and ensure it *does* use nuclear targets
		nuSQUIDS nus(marray<double, 1>({3},{100*units.GeV, 1e3*units.GeV, 1e4*units.GeV}), 3, neutrino, true);
		nus.SetNeutrinoCrossSections(wcg24Nuclear);
		std::string prem_path=getResourcePath()+"/astro/EARTH_MODEL_PREM_wIso.dat";
        auto earth=std::make_shared<EarthAtm>(prem_path);
        nus.Set_Body(earth);
        auto track=earth->MakeTrackWithCosine(-1);
        nus.Set_Track(std::make_shared<EarthAtm::Track>(track));
        nus.InitializeInteractions();
        auto intData=nus.GetInteractionStructure();
        assert(intData);
        assert(std::any_of(intData->targets.begin(), intData->targets.end(), isNuclearPDGCode));
	}
}
