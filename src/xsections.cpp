 /******************************************************************************
 *    This program is free software: you can redistribute it and/or modify     *
 *   it under the terms of the GNU General Public License as published by      *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   This program is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                             *
 *   Authors:                                                                  *
 *      Carlos Arguelles (University of Wisconsin Madison)                     *
 *         carguelles@icecube.wisc.edu                                         *
 *      Jordi Salvado (University of Wisconsin Madison)                        *
 *         jsalvado@icecube.wisc.edu                                           *
 *      Christopher Weaver (University of Wisconsin Madison)                   *
 *         chris.weaver@icecube.wisc.edu                                       *
 ******************************************************************************/

#include <nuSQuIDS/xsections.h>
#include <fstream>
#include <iostream>

#include "nuSQuIDS/AdaptiveQuad.h"
#include "nuSQuIDS/resources.h"

namespace{
	//Most of the time the cross section functions are _very_ smooth, and indeed
	//implemented as piece-wise cubic functions. This function therefore attempts 
	//to optimistically do an integral with a bare minimum of evaluations, 
	//starting with just two. If those values match to within the tolerance 
	//(treating them as 'quadrature rules' with degree of exactness zero), it  
	//treats the integrand as constant, otherwise it adds more evaluation points  
	//and uses quadrature rules which are exact for higher degress, namely one 
	//and three, returning the results of those rules if they agree to within 
	//the tolerance and otherwise continuing to escalate the sophistication of 
	//the approximation until finally calling the fully generaly adaptive
	//integration procedure. 
	//This function should not but used on functions which are substantially 
	//non-monotonic as it will, for example, return zero if asked to integrate
	//sin(x) over the interval [0,pi], regardless of tolerance. 
	template<typename FuncType>
	double fastInt(FuncType integrand, double a, double b, double tol){
		double abcissas[5]={-1,-1/sqrt(3),0,1/sqrt(3),1};
		const double midpoint=AdaptiveQuad::midpoint(a,b);
		const double halfWidth=(b-a)/2;
		auto checkAgreement=[](double est1, double est2, double tol){
			double a1=std::abs(est1);
			double a2=std::abs(est2);
			if(a1<a2) std::swap(a1,a2);
			return((a2 && tol>(a1/a2)-1) || (!a2 && a1<tol));
		};
		auto x=[&](unsigned int i){ return midpoint+abcissas[i]*halfWidth; };
		double y[5]; //array of integrand evaluations
		y[0]=integrand(x(0));
		y[4]=integrand(x(4));
		//check whether the integrand appers to be constant. This is somewhat
		//dangerous, as it can be fooled by the introduction of any term with 
		//even degree.  
		if(checkAgreement(y[0],y[4],tol))
			return halfWidth*(y[0]+y[4]); //actually exact for degree 1
		
		//use real quadrature rules of with degrees of exactness 3 and 1, 
		//allowing treatment of linear integrands. 
		y[1]=integrand(x(1));
		y[3]=integrand(x(3));
		double Legendre2=y[1]+y[3];
		double Lobatto2=y[0]+y[4]; //a.k.a. trapezoid rule
		if(checkAgreement(Legendre2,Lobatto2,tol))
			return halfWidth*Legendre2; //use result which is good to degree 3
		
		//try upgrading our calculation to be entirely exact for degree 3
		y[2]=integrand(x(2));
		double Lobatto3=(y[0]+4*y[2]+y[4])/3;
		if(checkAgreement(Legendre2,Lobatto3,tol))
			return halfWidth*Lobatto3;
		
		//otherwise, call the general adaptive integration, but give it our
		//known endpoint values to save recalculating them
		AdaptiveQuad::Options opt;
		opt.fa=y[0];
		opt.fb=y[4];
		return AdaptiveQuad::integrate(integrand,a,b,tol,&opt);
	}
}

namespace nusquids{

double NeutrinoCrossSections::AverageTotalCrossSection(double EMin, double EMax, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const{
	auto integrand=[=](double e_in)->double{
		return this->TotalCrossSection(e_in,flavor,neutype,current);
	};
	return fastInt(integrand,EMin,EMax,1e-3)/(EMax-EMin);
}
	
double NeutrinoCrossSections::AverageSingleDifferentialCrossSection(double E1, double E2Min, double E2Max, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const{
	auto integrand=[=](double e_out)->double{
		return this->SingleDifferentialCrossSection(E1,e_out,flavor,neutype,current);
	};
	return fastInt(integrand,E2Min,E2Max,1e-3)/(E2Max-E2Min);
}

std::tuple<AkimaSpline,double,double> NeutrinoDISCrossSectionsFromTables::read1DInterpolationFromText(const std::string& path){
	marray<double,2> rawData=quickread(path);
	if(rawData.extent(1)!=2)
		throw std::runtime_error(path+" does not contain two columns of data");
	const std::size_t nEntries=rawData.extent(0);
	std::vector<double> en, sigma;
	en.reserve(nEntries);
	sigma.reserve(nEntries);
	for(std::size_t i=0; i<nEntries; i++){
		en.push_back(log10(rawData[i][0]*GeV)); //note conversion to eV
		sigma.push_back(log10(rawData[i][1]));
	}
	double emin=rawData[0][0]*GeV;
	double emax=rawData[nEntries-1][0]*GeV;
	return std::make_tuple(AkimaSpline(en,sigma),emin,emax);
}

std::tuple<BiCubicInterpolator,double,double> NeutrinoDISCrossSectionsFromTables::read2DInterpolationFromText(const std::string& path){
	marray<double,2> rawData=quickread(path);
	if(rawData.extent(1)!=3)
		throw std::runtime_error(path+" does not contain three columns of data");
	double prevEn=0;
	std::size_t zSteps=0;
	for(; zSteps<rawData.extent(0); zSteps++){
		if(zSteps && rawData[zSteps][0]!=prevEn)
			break;
		else
			prevEn=rawData[zSteps][0];
	}
	std::size_t enSteps=rawData.extent(0)/zSteps;
	if(zSteps==rawData.extent(0) || enSteps*zSteps!=rawData.extent(0))
		throw std::runtime_error(path+" does not appear to represent a correctly ordered 2D table");
	
	marray<double,2> data({zSteps,enSteps});
	marray<double,1> xcoords({enSteps}), ycoords({zSteps});
	
	for(std::size_t j=0; j<enSteps; j++)
		xcoords[j]=log10(rawData[j*zSteps][0]*GeV); // note conversion to eV
	
	//TODO: here we throw away the actual data for the z values, and 
	//recompute them assuming they are linear. This assumption could be relaxed.
	for(std::size_t j=0; j<zSteps; j++)
		ycoords[j]=j/double(zSteps-1);
	
	//Constant for a small value in log-space, used to replace log(0). Making
	//this too small (or leaving infinities) causes the interpolation to misbehave. 
	const double tiny=-50;
	for(std::size_t j=0; j<enSteps; j++){
		for(std::size_t i=0; i<zSteps; i++){
			data[i][j]=log10(rawData[j*zSteps+i][2]);
			if(!std::isfinite(data[i][j]))
				data[i][j]=tiny;
		}
	}
	
	return std::make_tuple(BiCubicInterpolator(std::move(data),std::move(xcoords),std::move(ycoords)),
						   rawData[0][0]*GeV,rawData[rawData.extent(0)-1][0]*GeV);
}

bool NeutrinoDISCrossSectionsFromTables::isHDF(const std::string& path){
	std::ifstream infile(path);
	if(!infile)
		return false;
	unsigned char magic[8];
	infile.read((char*)magic,8);
	if(!infile)
		return false;
	return(magic[0]==0x89 && magic[1]=='H' && magic[2]=='D' && magic[3]=='F'
	  && magic[4]=='\r' && magic[5]=='\n' && magic[6]==0x1A && magic[7]=='\n');
}

NeutrinoDISCrossSectionsFromTables::NeutrinoDISCrossSectionsFromTables():
NeutrinoDISCrossSectionsFromTables(getResourcePath()+"/xsections/csms_square.h5"){}

NeutrinoDISCrossSectionsFromTables::NeutrinoDISCrossSectionsFromTables(std::string pathOrPrefix){
	if(isHDF(pathOrPrefix))
		ReadHDF(pathOrPrefix);
	else
		ReadText(pathOrPrefix);
}

double NeutrinoDISCrossSectionsFromTables::TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const{
	// we assume that sterile neutrinos are truly sterile
	if (not (flavor == electron or flavor == muon or flavor == tau))
		return 0;
	if (Enu > Emax)
		throw std::runtime_error("NeutrinoCrossSections::TotalCrossSection: Only DIS cross sections are included in a limited range. Interpolation re\
		quested below "+std::to_string(Emin/GeV)+" GeV or above "+std::to_string(Emax/GeV)+" GeV. E_nu = " + std::to_string(Enu/GeV) + " [GeV].");

	if (Enu < 10*GeV and Emin > 10*GeV) {
		std::clog << "NeutrinoCrossSections::TotalCrossSection: Neglecting the neutrino cross section below 10 GeV." << std::endl;
		return std::numeric_limits<double>::min();
	} else if (Enu < Emin) {
		// use approximate linear scaling below Emin GeV which is a good approximation for DIS up to 10 GeV
		return TotalCrossSection(Emin, flavor, neutype, current)*(Enu/Emin);
	}
	
	const AkimaSpline& sigma=
	(neutype==neutrino ? (current==CC ? s_CC_nu : s_NC_nu)
					   : (current==CC ? s_CC_nubar : s_NC_nubar));
	return pow(10.,sigma(log10(Enu)));
}

double NeutrinoDISCrossSectionsFromTables::SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const{
	// we assume that sterile neutrinos are truly sterile
	if (not (flavor == electron or flavor == muon or flavor == tau))
		return 0;

	if (E1 > Emax)
		throw std::runtime_error("NeutrinoCrossSections::SingleDifferentialCrossSection: Only DIS cross sections are included. Interpolation re\
		quested below "+std::to_string(Emin/GeV)+" GeV or above "+std::to_string(Emax/GeV)+" GeV. E_nu = " + std::to_string(E1/GeV) + " [GeV].");
	if (E1 <= Emin)
		return std::numeric_limits<double>::min();

	const BiCubicInterpolator& dsdy=
	(neutype==neutrino ? (current==CC ? dsdy_CC_nu : dsdy_NC_nu)
					   : (current==CC ? dsdy_CC_nubar : dsdy_NC_nubar));
	double z=(E2-Emin)/(E1-Emin);
	double val=dsdy(log10(E1),z);
	if(E2<=Emin){ //if extrapolating, make sure no blow-up happens
		double altVal=dsdy(log10(E1),0);
		if(val>altVal)
			val=altVal;
	}
	return pow(10.,val)/(E1/GeV);
}

double NeutrinoDISCrossSectionsFromTables::AverageSingleDifferentialCrossSection(double E1, double E2Min, double E2Max, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const{ 
	// we assume that sterile neutrinos are truly sterile
	if (not (flavor == electron or flavor == muon or flavor == tau))
		return 0;

	if (E1 > Emax)
		throw std::runtime_error("NeutrinoCrossSections::SingleDifferentialCrossSection: Only DIS cross sections are included. Interpolation re\
		quested below "+std::to_string(Emin/GeV)+" GeV or above "+std::to_string(Emax/GeV)+" GeV. E_nu = " + std::to_string(E1/GeV) + " [GeV].");
	if (E1 <= Emin)
		return std::numeric_limits<double>::min();
		
	if(E2Min>E2Max)
		std::swap(E2Min,E2Max);
	if(E2Min==E2Max)
		return SingleDifferentialCrossSection(E1, E2Max, flavor, neutype, current);
	
	auto integrand=[=](double e_out)->double{
		return this->SingleDifferentialCrossSection(E1,e_out,flavor,neutype,current);
	};
	return fastInt(integrand,E2Min,E2Max,1e-3)/(E2Max-E2Min);
}

void NeutrinoDISCrossSectionsFromTables::ReadText(const std::string& prefix){
	double eMinTmp, eMaxTmp;
	std::tie(s_CC_nu,eMinTmp,eMaxTmp)=read1DInterpolationFromText(prefix+"nu_sigma_CC.dat");
	Emin=eMinTmp;
	Emax=eMaxTmp;
	std::tie(s_NC_nu,eMinTmp,eMaxTmp)=read1DInterpolationFromText(prefix+"nu_sigma_NC.dat");
	if(eMinTmp!=Emin || eMaxTmp!=Emax)
		throw std::runtime_error(prefix+"nu_sigma_NC.dat has different energy domain than "+prefix+"nu_sigma_CC.dat");
	std::tie(s_CC_nubar,eMinTmp,eMaxTmp)=read1DInterpolationFromText(prefix+"nubar_sigma_CC.dat");
	if(eMinTmp!=Emin || eMaxTmp!=Emax)
		throw std::runtime_error(prefix+"nubar_sigma_CC.dat has different energy domain than "+prefix+"nu_sigma_CC.dat");
	std::tie(s_NC_nubar,eMinTmp,eMaxTmp)=read1DInterpolationFromText(prefix+"nubar_sigma_NC.dat");
	if(eMinTmp!=Emin || eMaxTmp!=Emax)
		throw std::runtime_error(prefix+"nubar_sigma_NC.dat has different energy domain than "+prefix+"nu_sigma_CC.dat");
	
	std::tie(dsdy_CC_nu,eMinTmp,eMaxTmp)=read2DInterpolationFromText(prefix+"nu_dsde_CC.dat");
	if(eMinTmp!=Emin || eMaxTmp!=Emax)
		throw std::runtime_error(prefix+"nu_dsde_CC.dat has different energy domain than "+prefix+"nu_sigma_CC.dat");
	std::tie(dsdy_NC_nu,eMinTmp,eMaxTmp)=read2DInterpolationFromText(prefix+"nu_dsde_NC.dat");
	if(eMinTmp!=Emin || eMaxTmp!=Emax)
		throw std::runtime_error(prefix+"nu_dsde_NC.dat has different energy domain than "+prefix+"nu_sigma_CC.dat");
	std::tie(dsdy_CC_nubar,eMinTmp,eMaxTmp)=read2DInterpolationFromText(prefix+"nubar_dsde_CC.dat");
	if(eMinTmp!=Emin || eMaxTmp!=Emax)
		throw std::runtime_error(prefix+"nubar_dsde_CC.dat has different energy domain than "+prefix+"nu_sigma_CC.dat");
	std::tie(dsdy_NC_nubar,eMinTmp,eMaxTmp)=read2DInterpolationFromText(prefix+"nubar_dsde_NC.dat");
	if(eMinTmp!=Emin || eMaxTmp!=Emax)
		throw std::runtime_error(prefix+"nubar_dsde_NC.dat has different energy domain than "+prefix+"nu_sigma_CC.dat");
}

void NeutrinoDISCrossSectionsFromTables::WriteText(const std::string& prefix) const{
	auto write1D=[this](const AkimaSpline& data, const std::string& path){
		std::ofstream outfile(path);
		if(!outfile)
			throw std::runtime_error("Unable to open "+path+" for writing");
		//this makes a copy, which is unfortunate, but not a disaster for the 1D data
		auto ab=data.getAbscissas();
		auto ord=data.getOrdinates();
		for(std::size_t i=0, n=ab.extent(0); i!=n; i++)
			outfile << pow(10.,ab[i])/GeV << ' ' << pow(10.,ord[i]) << '\n';
	};
	auto write2D=[this](const BiCubicInterpolator& data, const std::string& path){
		std::ofstream outfile(path);
		if(!outfile)
			throw std::runtime_error("Unable to open "+path+" for writing");
		const auto& x=data.getXCoords();
		const auto& y=data.getYCoords();
		const auto& z=data.getData();
		for(std::size_t i=0, iMax=x.extent(0); i!=iMax; i++){
			double en=pow(10.,x[i])/GeV;
			for(std::size_t j=0, jMax=y.extent(0); j!=jMax; j++)
				outfile << en << ' ' << y[j] << ' ' << pow(10.,z[i][j]) << '\n';
		}
	};
	write1D(s_CC_nu,prefix+"nu_sigma_CC.dat");
	write1D(s_NC_nu,prefix+"nu_sigma_NC.dat");
	write1D(s_CC_nubar,prefix+"nubar_sigma_CC.dat");
	write1D(s_NC_nubar,prefix+"nubar_sigma_NC.dat");
	write2D(dsdy_CC_nu,prefix+"nu_dsde_CC.dat");
	write2D(dsdy_NC_nu,prefix+"nu_dsde_NC.dat");
	write2D(dsdy_CC_nubar,prefix+"nubar_dsde_CC.dat");
	write2D(dsdy_NC_nubar,prefix+"nubar_dsde_NC.dat");
}

void NeutrinoDISCrossSectionsFromTables::ReadHDF(const std::string& path){
	try{
	H5File h5file(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
	
	marray<double,1> energies,zs;
	readArrayH5(h5file, "energies", energies);
	readArrayH5(h5file, "zs", zs);
	Emin=pow(10.,energies.front());
	Emax=pow(10.,energies.back());
	
	marray<double,1> buffer1;
	auto readTotal=[&](AkimaSpline& dest, const std::string& tableName){
		readArrayH5(h5file, tableName, buffer1);
		if(buffer1.extent(0)!=energies.extent(0))
			throw std::runtime_error(tableName+" does not have the same dimensions as energies");
		dest=AkimaSpline(energies.get_data(),buffer1.get_data(),energies.size());
	};
	readTotal(s_CC_nu,"s_CC_nu");
	readTotal(s_NC_nu,"s_NC_nu");
	readTotal(s_CC_nubar,"s_CC_nubar");
	readTotal(s_NC_nubar,"s_NC_nubar");
	
	marray<double,2> buffer2;
	auto readDifferential=[&](BiCubicInterpolator& dest, const std::string& tableName){
		readArrayH5(h5file, tableName, buffer2); //reallocates the target array 
		if(buffer2.extent(0)!=energies.extent(0))
			throw std::runtime_error(tableName+" does not have the same first dimension as energies");
		if(buffer2.extent(1)!=zs.extent(0))
			throw std::runtime_error(tableName+"'s second dimension does not have the same size as zs");
		dest=BiCubicInterpolator(std::move(buffer2),energies,zs);
	};
	readDifferential(dsdy_CC_nu,"dsdy_CC_nu");
	readDifferential(dsdy_NC_nu,"dsdy_NC_nu");
	readDifferential(dsdy_CC_nubar,"dsdy_CC_nubar");
	readDifferential(dsdy_NC_nubar,"dsdy_NC_nubar");
	}
	catch(std::runtime_error& err){
		throw std::runtime_error("Failed to read cross sections from "+path+": "+err.what());
	}
}

void NeutrinoDISCrossSectionsFromTables::WriteHDF(const std::string& path, unsigned int compressionLevel) const{
	using StoreType=float;
	H5File h5file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
	writeArrayH5<StoreType>(h5file, "energies", dsdy_CC_nu.getXCoords(), compressionLevel);
	writeArrayH5<StoreType>(h5file, "zs", dsdy_CC_nu.getYCoords(), compressionLevel);
	writeArrayH5<StoreType>(h5file, "s_CC_nu", s_CC_nu.getOrdinates(), compressionLevel);
	writeArrayH5<StoreType>(h5file, "s_NC_nu", s_NC_nu.getOrdinates(), compressionLevel);
	writeArrayH5<StoreType>(h5file, "s_CC_nubar", s_CC_nubar.getOrdinates(), compressionLevel);
	writeArrayH5<StoreType>(h5file, "s_NC_nubar", s_NC_nubar.getOrdinates(), compressionLevel);
	writeArrayH5<StoreType>(h5file, "dsdy_CC_nu", dsdy_CC_nu.getData(), compressionLevel);
	writeArrayH5<StoreType>(h5file, "dsdy_NC_nu", dsdy_NC_nu.getData(), compressionLevel);
	writeArrayH5<StoreType>(h5file, "dsdy_CC_nubar", dsdy_CC_nubar.getData(), compressionLevel);
	writeArrayH5<StoreType>(h5file, "dsdy_NC_nubar", dsdy_NC_nubar.getData(), compressionLevel);
}

double NeutrinoDISCrossSectionsFromTables_V1::LinInter(double x,double xM, double xP,double yM,double yP) const{
  return yM + (yP-yM)*(x-xM)/(xP-xM);
}

///do linear interpolation on a triangle
double TriLinInter(double x, double y, 
                   double x1, double y1, double z1,
                   double x2, double y2, double z2,
                   double x3, double y3, double z3){
 	//convert to barycentric coordinates
 	double denom=(y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
 	double l1=((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/denom;
 	double l2=((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/denom;
 	double l3=1-(l1+l2);
 	return l1*z1 + l2*z2 + l3*z3;
}

double NeutrinoDISCrossSectionsFromTables_V1::TotalCrossSection(double Enu, NeutrinoFlavor flavor,
                           NeutrinoType neutype, Current current) const{
  // we assume that sterile neutrinos are truly sterile
  if (not (flavor == electron or flavor == muon or flavor == tau))
    return 0.0;

  if (Enu > Emax)
    throw std::runtime_error("NeutrinoCrossSections::TotalCrossSection: Only DIS cross sections are included in a limited range. Interpolation re\
quested below "+std::to_string(Emin/GeV)+" GeV or above "+std::to_string(Emax/GeV)+" GeV. E_nu = " + std::to_string(Enu/GeV) + " [GeV].");

  if (Enu < 10*GeV and Emin > 10*GeV) {
    std::clog << "NeutrinoCrossSections::TotalCrossSection: Neglecting the neutrino cross section below 10 GeV." << std::endl;
    return std::numeric_limits<double>::min();
  } else if (Enu < Emin) {
    // use approximate linear scaling below Emin GeV which is a good approximation for DIS up to 10 GeV
    return TotalCrossSection(Emin, flavor, neutype, current)*(Enu/Emin);
  }

  // convert to GeV
  Enu /= GeV;

  double logE = log(Enu);
  double dlogE = logE_data_range[1]-logE_data_range[0];
  size_t idx = static_cast<size_t>((logE-logE_data_range[0])/dlogE);
  if(idx==div)
    idx--;
  const marray<double,3>& sigma=(current==CC ? s_CC_data : s_NC_data);
  return(LinInter(logE,logE_data_range[idx],logE_data_range[idx+1],sigma[neutype][flavor][idx],sigma[neutype][flavor][idx+1]));
}

double NeutrinoDISCrossSectionsFromTables_V1::SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const{
  // we assume that sterile neutrinos are trully sterile
  if (not (flavor == electron or flavor == muon or flavor == tau))
    return 0.0;

  if (E1 > Emax)
    throw std::runtime_error("NeutrinoCrossSections::SingleDifferentialCrossSection: Only DIS cross sections are included. Interpolation re\
quested below "+std::to_string(Emin/GeV)+" GeV or above "+std::to_string(Emax/GeV)+" GeV. E_nu = " + std::to_string(E1/GeV) + " [GeV].");
  if (E1 < Emin)
    return std::numeric_limits<double>::min();

  // convert to GeV
  E1 /= GeV;
  E2 /= GeV;

  double logE1 = log(E1);
  double logE2 = log(E2);
  double dlogE = logE_data_range[1]-logE_data_range[0];

  size_t loge_M1 = static_cast<size_t>((logE1-logE_data_range[0])/dlogE);
  size_t loge_M2 = static_cast<size_t>((logE2-logE_data_range[0])/dlogE);

  if ((loge_M2 > div-1) or (loge_M1 > div-1) or (E2 >= E1))
    return 0.0;

  //std::cout << E1 << " " << E2 << " " << loge_M1 << " " << loge_M2 << " " << div << std::endl;
  double phiMM,phiMP,phiPM,phiPP;
  if (current == CC){
    phiMM = dsde_CC_data[neutype][flavor][loge_M1][loge_M2];
    phiMP = dsde_CC_data[neutype][flavor][loge_M1][loge_M2+1];
    if ( loge_M1 == div-1 ){
      // we are at the boundary, cannot bilinearly interpolate
      return LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],
          phiMM,phiMP);
    }
    phiPM = dsde_CC_data[neutype][flavor][loge_M1+1][loge_M2];
    phiPP = dsde_CC_data[neutype][flavor][loge_M1+1][loge_M2+1];
  } else if (current == NC){
    phiMM = dsde_NC_data[neutype][flavor][loge_M1][loge_M2];
    phiMP = dsde_NC_data[neutype][flavor][loge_M1][loge_M2+1];
    if ( loge_M1 == div-1 ){
      // we are at the boundary, cannot bilinearly interpolate
      return LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],
          phiMM,phiMP);
    }
    phiPM = dsde_NC_data[neutype][flavor][loge_M1+1][loge_M2];
    phiPP = dsde_NC_data[neutype][flavor][loge_M1+1][loge_M2+1];
  } else
    throw std::runtime_error("nuSQUIDS::XSECTIONS::ERROR::Current type unkwown.");

  //If close to the E2==E1 boundary, omit the non-physical E2>E1 point from the interpolation
  if(loge_M1==loge_M2 && phiMP==0){
    return TriLinInter(logE1, logE2, 
                       logE_data_range[loge_M1], logE_data_range[loge_M2], phiMM,
                       logE_data_range[loge_M1+1], logE_data_range[loge_M2], phiPM,
                       logE_data_range[loge_M1+1], logE_data_range[loge_M2+1], phiPP);
  }
    
  return LinInter(logE1,logE_data_range[loge_M1],logE_data_range[loge_M1+1],
           LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],phiMM,phiMP),
           LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],phiPM,phiPP));
}

void NeutrinoDISCrossSectionsFromTables_V1::ReadText(std::string root){
       std::string filename_dsde_CC = root+"dsde_CC.dat";
       std::string filename_dsde_NC = root+"dsde_NC.dat";
       std::string filename_sigma_CC = root+"sigma_CC.dat";
       std::string filename_sigma_NC = root+"sigma_NC.dat";

       // check if files exist for this energies and divisions
       if(
          fexists(filename_dsde_CC) and
          fexists(filename_dsde_NC) and
          fexists(filename_sigma_CC) and
          fexists(filename_sigma_NC)
          )
       {
          // read data tables
          marray<double,2> dsde_CC_raw_data = quickread(filename_dsde_CC);
          marray<double,2> dsde_NC_raw_data = quickread(filename_dsde_NC);
          marray<double,2> sigma_CC_raw_data = quickread(filename_sigma_CC);
          marray<double,2> sigma_NC_raw_data = quickread(filename_sigma_NC);

          // check table shapes and get the number of energy nodes
          unsigned int data_e_size = 0;
          if( sigma_CC_raw_data.extent(0) == sigma_NC_raw_data.extent(0) and
              sigma_NC_raw_data.extent(1) == sigma_NC_raw_data.extent(1) )
            data_e_size = sigma_CC_raw_data.extent(0);
          else
            throw std::runtime_error("nuSQUIDS::xsections::init: Data tables not the same size.");

          // getting the raw data energy node values
          logE_data_range.resize(data_e_size);
          for( int ie = 0; ie < data_e_size; ie ++)
            logE_data_range[ie] = log(sigma_CC_raw_data[ie][0]);

          Emin = sigma_CC_raw_data[0][0]*GeV;
          Emax = sigma_CC_raw_data[data_e_size-1][0]*GeV;
          div = data_e_size;

          // convert raw data tables into formatted marrays
          s_CC_data.resize(std::vector<size_t>{2,3,data_e_size});
          s_NC_data.resize(std::vector<size_t>{2,3,data_e_size});
          dsde_CC_data.resize(std::vector<size_t>{2,3,data_e_size,data_e_size});
          dsde_NC_data.resize(std::vector<size_t>{2,3,data_e_size,data_e_size});
          for(NeutrinoType neutype : {neutrino,antineutrino}){
            for(NeutrinoFlavor flavor : {electron,muon,tau}){
              for(unsigned int e1 = 0; e1 < data_e_size; e1++){
                s_CC_data[neutype][flavor][e1] = sigma_CC_raw_data[e1][1+2*((int)flavor)+(int)neutype];
                s_NC_data[neutype][flavor][e1] = sigma_NC_raw_data[e1][1+2*((int)flavor)+(int)neutype];
                for (unsigned int e2 = 0; e2 < data_e_size; e2++){
                  dsde_CC_data[neutype][flavor][e1][e2] = dsde_CC_raw_data[e1*data_e_size+e2][2+2*(static_cast<int>(flavor))+static_cast<int>(neutype)];
                  dsde_NC_data[neutype][flavor][e1][e2] = dsde_NC_raw_data[e1*data_e_size+e2][2+2*(static_cast<int>(flavor))+static_cast<int>(neutype)];
                }
              }
            }
          }
  } else {
    throw std::runtime_error("nuSQUIDS::XSECTIONS::ERROR::Cross section files not found.");
  }
  // declare the object as initialized;
  is_init = true;
}
  
NeutrinoDISCrossSectionsFromTables_V1::NeutrinoDISCrossSectionsFromTables_V1():
  NeutrinoDISCrossSectionsFromTables_V1(getResourcePath()+"/xsections/csms.h5"){
}
    
NeutrinoDISCrossSectionsFromTables_V1::NeutrinoDISCrossSectionsFromTables_V1(std::string path){
  //If a single file, read HDF5
  if(fexists(path)){
    {
      H5File h5file(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      readH5Attribute(h5file, "Emin", Emin);
      readH5Attribute(h5file, "Emax", Emax);
      readArrayH5(h5file, "s_CC", s_CC_data);
      readArrayH5(h5file, "s_NC", s_NC_data);
      readArrayH5(h5file, "dsDE_CC", dsde_CC_data);
      readArrayH5(h5file, "dsDE_NC", dsde_NC_data);
    }
    //TODO: make sanity checking more user friendly
    assert(dsde_CC_data.extent(2)==dsde_CC_data.extent(3));
    assert(dsde_NC_data.extent(2)==dsde_NC_data.extent(3));
    assert(dsde_CC_data.extent(2)==dsde_NC_data.extent(2));
    assert(dsde_CC_data.extent(0)==2);
    assert(dsde_NC_data.extent(0)==2);
    assert(dsde_CC_data.extent(1)==3);
    assert(dsde_NC_data.extent(1)==3);
    
    div=dsde_CC_data.extent(2);
    unsigned int data_e_size=div;
    auto raw_logE_data_range=logspace(Emin/GeV,Emax/GeV,div);
    logE_data_range.resize(data_e_size);
    std::transform(raw_logE_data_range.begin(),raw_logE_data_range.end(),logE_data_range.begin(),
                   (double(*)(double))log);
    is_init=true;
  }
  else //otherwise, try to read several text files
    ReadText(path);
}

NeutrinoDISCrossSectionsFromTables_V1::~NeutrinoDISCrossSectionsFromTables_V1(){}
    
void NeutrinoDISCrossSectionsFromTables_V1::WriteHDF(std::string path) const{
  H5File h5file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  addH5Attribute(h5file, "Emin", Emin);
  addH5Attribute(h5file, "Emax", Emax);
  writeArrayH5<double>(h5file, "s_CC", s_CC_data, 0);
  writeArrayH5<double>(h5file, "s_NC", s_NC_data, 0);
  writeArrayH5<double>(h5file, "dsDE_CC", dsde_CC_data, 0);
  writeArrayH5<double>(h5file, "dsDE_NC", dsde_NC_data, 0);
}
  
void NeutrinoDISCrossSectionsFromTables_V1::WriteText(std::string basePath) const{
  std::string filename_sigma_CC = basePath+"sigma_CC.dat";
  std::string filename_sigma_NC = basePath+"sigma_NC.dat";
  std::string filename_dsde_CC = basePath+"dsde_CC.dat";
  std::string filename_dsde_NC = basePath+"dsde_NC.dat";
  
  auto writeTotal=[](const std::string& path, const marray<double,3>& data,
                     const std::vector<double>& logEnergies, double GeV){
    std::ofstream out(path);
    for(size_t ie=0; ie<logEnergies.size(); ie++){
      out << exp(logEnergies[ie]) << ' ';
      for(size_t fl : {0,1,2}){
        for(size_t pt : {0,1})
          out << data[pt][fl][ie] << ' ';
      }
      out << '\n';
    }
  };
  writeTotal(filename_sigma_CC,s_CC_data,logE_data_range,GeV);
  writeTotal(filename_sigma_NC,s_NC_data,logE_data_range,GeV);
  
  auto writeDifferential=[](const std::string& path, const marray<double,4>& data,
                            const std::vector<double>& logEnergies, double GeV){
    std::ofstream out(path);
    for(size_t ie1=0; ie1<logEnergies.size(); ie1++){
      double e1=exp(logEnergies[ie1]);
      for(size_t ie2=0; ie2<logEnergies.size(); ie2++){
        double e2=exp(logEnergies[ie2]);
        out << e1 << ' ' << e2 << ' ';
        for(size_t fl : {0,1,2}){
          for(size_t pt : {0,1})
            out << data[pt][fl][ie1][ie2] << ' ';
        }
        out << '\n';
      }
    }
  };
  writeDifferential(filename_dsde_CC,dsde_CC_data,logE_data_range,GeV);
  writeDifferential(filename_dsde_NC,dsde_NC_data,logE_data_range,GeV);
}
  
GlashowResonanceCrossSection::GlashowResonanceCrossSection():
fermi_scale(pow(constants.GF/constants.cm, 2)*constants.electron_mass/constants.pi),
M_W(80.385*constants.GeV),
W_total(2.085*constants.GeV)
{}
  
//with no non-constant data, copy construction is the same as default construction
GlashowResonanceCrossSection::GlashowResonanceCrossSection(const GlashowResonanceCrossSection&):
GlashowResonanceCrossSection(){}
//same for move construction
GlashowResonanceCrossSection::GlashowResonanceCrossSection(GlashowResonanceCrossSection&&):
GlashowResonanceCrossSection(){}
  
GlashowResonanceCrossSection::~GlashowResonanceCrossSection(){}
  
double GlashowResonanceCrossSection::TotalCrossSection(double Enu, NeutrinoFlavor flavor,
 NeutrinoType neutype, Current current) const{
  // only treat the glashow resonance
  if (not (flavor == electron and neutype == antineutrino and current == GR))
    return 0;
  // calculate total cross-section for nuebar + e- -> numubar + mu-, then divide
  // by the W- -> numubar + mu- branching fraction to get total cross-section
  // see e.g. arxiv:1108.3163v2, Eq. 2.1
  const double &m = constants.electron_mass;
  const double &mu = constants.muon_mass;
  return 2*fermi_scale*Enu*pow(1. - (mu*mu - m*m)/(2*m*Enu), 2)
   / (pow(1. - 2*m*Enu/(M_W*M_W), 2) + pow(W_total/M_W, 2)) / B_Muon / 3;
}
  
double GlashowResonanceCrossSection::SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const{
  // only treat the glashow resonance
  if (not (flavor == electron and neutype == antineutrino and current == GR))
    return 0;
  if (E2 > E1)
    return 0;
  // differential cross section for leptonic final states only
  // NB: assumes that the branching fractions are identical for all 3 families
  double xl = E2/E1;
  return B_Muon*3*TotalCrossSection(E1, flavor, neutype, current)*(xl*xl)/E1*constants.GeV;
}
  
double GlashowResonanceCrossSection::DoubleDifferentialCrossSection(double E, double x, double y, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const{
  return 0;
}
  
//K. Olive et al. (PDG), Chin. Phys. C38, 090001 (2014)
double GlashowResonanceCrossSection::B_Electron = .1071; //unc. .0016
double GlashowResonanceCrossSection::B_Muon     = .1063; //unc. .0015
double GlashowResonanceCrossSection::B_Tau      = .1138; //unc. .0021
double GlashowResonanceCrossSection::B_Hadronic = .6741; //unc. .0027

std::shared_ptr<const NeutrinoCrossSections> CrossSectionLibrary::crossSectionForTarget(PDGCode target) const{
    auto it = data.find(target);
    if(it==data.end())
        return {};
    return it->second;
}
    
bool CrossSectionLibrary::hasTarget(PDGCode target) const{
    auto it = data.find(target);
    return(it!=data.end());
}
    
CrossSectionLibrary loadDefaultCrossSections(){
    CrossSectionLibrary lib;
    
    std::string xsdir = getResourcePath()+"/xsections/";
    //old, isoscalar table
    //lib.addTarget(isoscalar_nucleon, NeutrinoDISCrossSectionsFromTables(XSECTION_LOCATION "csms_square.h5"));
    //shiny, new, per-target tables
    lib.addTarget(proton, NeutrinoDISCrossSectionsFromTables(xsdir+"csms_proton.h5"));
    lib.addTarget(neutron,NeutrinoDISCrossSectionsFromTables(xsdir+"csms_neutron.h5"));
    
    lib.addTarget(electron,GlashowResonanceCrossSection());
    return lib;
}
    
} // close namespace
