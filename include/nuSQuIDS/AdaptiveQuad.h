/*  Copyright 2020 C. Weaver
 
 Redistribution and use in source and binary forms, with or without modification, 
 are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, 
 this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, 
 this list of conditions and the following disclaimer in the documentation 
 and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace AdaptiveQuad{
	
	///\pre a<=b && isfinite(a) && isfinite(b)
	template<typename N>
	N midpoint(N a, N b){
		//Following dx.doi.org/10.1145/2493882 relaxing the precondition would 
		//require doing the following, but it is expensive and rarely useful:
		//if(a>b)
		//	std::swap(a,b);
		//if(a==-b)
		//	return N(0);
		//if(a==-std::numeric_limits<N>::infinity())
		//	return -std::numeric_limits<N>::max();
		//if(b==std::numeric_limits<N>::infinity())
		//	return std::numeric_limits<N>::max();
		return (a-a/2)+b/2;
	}
	
	struct Interval{
		///left endpoint
		double a;
		///right endpoint
		double b;
		///integrand evaluated at left endpoint
		double fa;
		///integrand evaluated at right endpoint
		double fb;
		///left split point
		double spa;
		///right split point
		double spb;
		///integrand evaluated at left split point
		double fspa;
		///integrand evaluated at right split point
		double fspb;
		///integral estimate
		double value;
		///second integral estimate
		double i2;
		
		Interval(double a, double b, double fa, double fb):
		a(a),b(b),fa(fa),fb(fb),value(std::numeric_limits<double>::quiet_NaN()){}
		
		template <typename Func>
		double evaluate(Func&& f, const bool full=false){
			const static double alpha=sqrt(2./3);
			const static double beta=1/sqrt(5);
			const static double x1=.942882415695480;
			const static double x2=.641853342345781;
			const static double x3=.236383199662150;
			const static double GanderGautschiAbscissas[13]={-1,-x1,-alpha,-x2,
			                                                 -beta,-x3,0,x3,beta,
			                                                 x2,alpha,x1,1};
			
			const double midpoint=AdaptiveQuad::midpoint(a,b);
			const double halfWidth=(b-a)/2;
			//function to scale the generic abscissas to this particular domain
			auto x=[&](unsigned int i)->double{ 
				return midpoint+GanderGautschiAbscissas[i]*halfWidth; 
			};
			
			double y[13]; //array of integrand evaluations
			//insert the endpoint values which are already known
			y[0]=fa;
			y[12]=fb;
			//evaluate interior points
			//odd numbered evaluates are only needed fo the 13 point extension
			for(unsigned int i=(full?1:2); i!=12; i+=(full?1:2))
				y[i]=f(x(i));
			//pull out the values which can be reused if this interval is subdivided
			spa=x(4);
			fspa=y[4];
			spb=x(8);
			fspb=y[8];
			//evaluate 4-point Gauss-Lobatto quadrature
			i2=(halfWidth/6.)*(y[0]+y[12]+5.*(y[4]+y[8]));
			//evaluate 7-point Kronrod extension
			value=(halfWidth/1470.)*(77.*(y[0]+y[12])+
			                         432.*(y[2]+y[10])+
			                         625.*(y[4]+y[8])+
			                         672.*y[6]);
			if(!full)
				return 0;
			//evaluate 13-point Kronrod extension
			return halfWidth*(.0158271919734802*(y[0]+y[12])+
			                  .0942738402188500*(y[1]+y[11])+
			                  .1550719873365850*(y[2]+y[10])+
			                  .1888215739601820*(y[3]+y[9])+
			                  .1997734052268590*(y[4]+y[8])+
			                  .224926465333340*(y[5]+y[7])+
			                  .242611071901408*y[6]);
		}
		
		bool canSubdivide() const{
			//Subdivision involves spliting into three subintervals, and 
			//evaluting five points in the interior of each, to this interval 
			//must be wide enough for all 1+5+1+5+1+5+1=19 points to be distinct.
			//Indeed, this is probably still not sufficient as we want these 
			//points to be unevenly distrbuted in a precise way, but the main
			//goal of this check is just to prevent runaway splitting. 
			//This serves the same purpose as the traditional 
			//`a<midpoint && midpoint<b` check. 
			//Here, 19 is rounded up to 32, largely for aesthetics. 
			return (b-a) > 
				32*std::numeric_limits<double>::epsilon()*std::max(std::abs(a),std::abs(b));
		}
		
		///Eavluate how bad this interval's uncertainty is relative to a tolerance
		double badness(double tol) const{
			return std::abs(value-i2)-tol;
		}
	};
	
	///This type represents optional inputs an outputs for integrate. 
	///fa and fb may be set before calling integrate if they happen to already 
	///be known, for example because multiple integrals over adjacent 
	///intervals are being computed. 
	///If passed an Options object, integrate will fill in fa and fb for 
	///possible reuse if they were not already supplied by the caller. 
	///It will set outOfTolerance if it is not able to meet its tolerance goal
	///by continued subdivision (if a subinterval with unacceptable uncertainty
	///becomes too small to further subdivide), and it will set uncertainty to
	///its estimate of the error on its result. 
	struct Options{
	public:
		Options():
		fa(std::numeric_limits<double>::quiet_NaN()),
		fb(std::numeric_limits<double>::quiet_NaN()),
		outOfTolerance(false)
		{}
	
		double fa,fb;
		bool inTolerance() const{ return !outOfTolerance; }

		bool outOfTolerance;
		double uncertainty;
	};
	
	///Integrate a function with adaptive gaussian quadrature. 
	///This is a modified version of the standard Gand and Gautschi algorithm
	///which splits intervals into three subintervals at a time rather than six.
	///\param integrand the function to be integrated
	///\param a the left boundary of the interval over which to integrate
	///\param b the right boundary of the interval over which to integrate
	///\param tol the relative error to be sought in computing the integral 
	///           estimate
	///\param options a pointer to an Options object which can convey additional
	///               inputs and outputs. See the Options type for details. 
	///\pre a and b should be finite
	template <typename Func>
	double integrate(Func&& integrand, double a, double b, double tol=1e-6, 
	                 Options* options=nullptr){
		if(tol<std::numeric_limits<double>::epsilon())
			throw std::runtime_error("Specified tolerance too small");
		
		//ensure that the interval is left to right, to simplify internals
		bool reverse=false;
		if(a>b){
			std::swap(a,b);
			reverse=true;
		}
		
		//compute function endpoint values if not given
		double fa,fb;
		if(options && !std::isnan(options->fa))
			fa=options->fa;
		else{
			fa=integrand(a);
			if(options) options->fa=fa;
		}
		if(options && !std::isnan(options->fb))
			fb=options->fb;
		else{
			fb=integrand(b);
			if(options) options->fb=fb;
		}
		
		//create initial interval with function endpoint values
		Interval i(a,b,fa,fb);
		if(!i.canSubdivide())
			throw std::runtime_error("Interval too narrow");
		
		double is,is1;
		is=is1=i.evaluate(integrand,true);
		double erri1=std::abs(i.value-is);
		double erri2=std::abs(i.i2-is);
		double R=(erri2!=0) ? erri1/erri2 : 1;
		if(R>0 && R<1)
			tol/=R;
		if(is==0)
			is=b-a;
		is=std::abs(is);
		tol*=is; //convert to an absolute tolerance
		if(i.badness(tol) <= 0){
			if(options){
				options->uncertainty=erri1;
				options->outOfTolerance=false;
			}
			//if we met our condition without subdividing, 
			//optimistically use the high order estimate, since we have it
			return (reverse?-is1:is1);
		}
		
		std::vector<Interval> intervals;
		intervals.reserve(10);
		intervals.push_back(i);
		bool outOfTolerance=false;
		double total=0, comp=0, unc=0;
		auto accumulate=[&total,&comp,&unc](double contribution, double est2){
			double y=contribution-comp;
			double t=total+y;
			comp=(t-total)-y;
			total=t;
			//Not being particularly concerned about the uncertainty on the 
			//uncertainty, we do not use compensated summation for it. 
			double uncc=contribution-est2;
			unc+=uncc*uncc;
		};
		while(!intervals.empty()){
			const Interval& i=intervals.back();
			Interval i1(i.a,i.spa,i.fa,i.fspa);
			Interval i2(i.spa,i.spb,i.fspa,i.fspb);
			Interval i3(i.spb,i.b,i.fspb,i.fb);
			intervals.pop_back(); //invalidates i
			i1.evaluate(integrand);
			i2.evaluate(integrand);
			i3.evaluate(integrand);
			//For each sub-interval, check whether there is any need to 
			//subdivide it further. If there is not, or subdivision is not 
			//possible, commit its contribution to the total estimate. 
			auto processSubinterval=[&](Interval& i){
				if(i.badness(tol)<=0 || !i.canSubdivide()){
					if(!i.canSubdivide())
						outOfTolerance=true;
					accumulate(i.value,i.i2);
				}
				else
					intervals.push_back(i);
			};
			processSubinterval(i1);
			processSubinterval(i2);
			processSubinterval(i3);
		}
		if(options){
			options->uncertainty=sqrt(unc)*(R>0 && R<1?R:1);
			options->outOfTolerance=outOfTolerance;
		}
		if(reverse)
			total=-total;
		return total;
	}
}
