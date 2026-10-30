#include <iostream>
#include <vector>

#include <SQuIDS/const.h>
#include <nuSQuIDS/body.h>

int main(){
	using namespace nusquids;
	squids::Const units;
	std::cout.precision(6);
	
	//define a rather peculiar Earth density profile which has simple but substantial changes very
	//close to the center and surface
	std::vector<double> fractionalRadii={0.0, 0.001, 0.002, 0.005, 0.01, 0.5, 0.99, 0.995, 0.998, 0.999, 1.0};
	std::vector<double> densities=      {1.0, 2.0,   3.0,   6.0,   11.0, 11.0,11.0, 6.0,   3.0,   2.0,   1.0};
	std::vector<double> efrac=          {0.0, 0.625, 0.125, 0.25,  0.5,  0.5, 0.5,  0.25,  0.125, 0.625, 0.0};
	
	EarthAtm earth(fractionalRadii, densities, efrac);
	//make a track which begins at the top of the atmosphere on the opposite side of the Earth,
	//passes through the full diameter, and stops when it reaches the surface again
	auto track=earth.MakeTrackWithCosine(-1);
	
	double atmHeight=earth.GetAtmosphereHeight()*units.km;
	double solidRadius=earth.GetRadius()*units.km;
	
	std::cout << "Track length: " << track.GetFinalX()-track.GetInitialX() << '\n';
	std::cout << "Expected track length: " << atmHeight+2*solidRadius << '\n';
	assert(track.GetFinalX()-track.GetInitialX() == atmHeight+2*solidRadius);
	
	//check that ten points in the atmosphere have small densities
	std::cout << "Atmosphere:\n";
	for(int i=0; i<10; i++){
		double x=i*atmHeight/10;
		track.SetX(x);
		double density=earth.density(track);
		std::cout << (atmHeight-x)/units.km << ' ' << density << std::endl;
		if(density>0.0012)
			std::cout << " Atmosphere density too large\n";
	}
	
	const double errTol=1e-6;
	
	//step through the first half of the path through the Earth, though in reverse (from the center
	//to the surface), ensuring that the evaluated density at each knot is equal to the input value
	std::cout << "In Earth, first half:\n";
	unsigned int idx=0;
	for(double u : fractionalRadii){
		double x=atmHeight+u*solidRadius;
		track.SetX(x);
		double density=earth.density(track);
		std::cout << u << ' ' << density << std::endl;
		double relErr=std::abs(density-densities[idx])/density;
		if(relErr>errTol)
			std::cout << " Too large relative error: " << relErr << std::endl;
		idx++;
	}
	//repeat for the second half of the Earth, again from center to surface, so no longer reversed
	std::cout << "In Earth, second half:\n";
	idx=0;
	for(double u : fractionalRadii){
		double x=atmHeight+solidRadius+u*solidRadius;
		track.SetX(x);
		double density=earth.density(track);
		std::cout << u << ' ' << density << std::endl;
		double relErr=std::abs(density-densities[idx])/density;
		if(relErr>errTol)
			std::cout << " Too large relative error: " << relErr << std::endl;
		idx++;
	}
}