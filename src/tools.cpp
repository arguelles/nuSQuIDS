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

#include <nuSQuIDS/tools.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

namespace nusquids{

bool fexists(std::string filepath)
{
  std::ifstream ifile(filepath.c_str());
  return static_cast<bool>(ifile);
}

marray<double,2> quickread(std::string filepath){
    // create and open file stream
    std::ifstream infile(filepath.c_str());

    if(!infile){
        throw std::runtime_error("Error: file could not be opened. Filepath " + filepath);
    }

    std::vector< std::vector<double> > table;
    std::string line;
    size_t column_number = -1;
    while(getline(infile,line)){
        std::vector<double> row;
        std::stringstream linestream(line);

        double data;
        while(linestream >> data){
            row.push_back(data);
        }
        if (!row.empty()){
            #ifdef quickread_DEBUG
            y = row.size();
            #endif
            if(!table.empty() and column_number != row.size())
              throw std::runtime_error("nuSQuIDS::tools::quickread: Number of columns in file not equal.");
            table.push_back(row);
        }

        column_number = row.size();
    }

    marray<double,2> otable {table.size(),column_number};
    for (unsigned int i = 0; i < table.size(); i++)
      for (unsigned int j = 0; j < column_number; j++)
        otable[i][j] = table[i][j];

    #ifdef quickread_DEBUG
    x = table.size();
    cout << "x: " << x << " y: "<< y << endl;
    cout << table[10][0] << endl;
    #endif
    return otable;
}

int quickwrite(std::string filepath, marray<double,2>& tbl){
    // create and open file stream
    std::ofstream outfile(filepath.c_str());

    if(!outfile)
        throw std::runtime_error("Error: file could not be created: " + filepath);

    outfile.precision(15);
    for (unsigned int i=0; i < tbl.extent(0); i++){
      for ( unsigned int j=0; j < tbl.extent(1); j++){
         outfile << tbl[i][j] << " ";
      }
      outfile << std::endl;
    }
    return 0;
}

marray<double,1> linspace(double Emin,double Emax,unsigned int div){
    if(div==0)
        throw std::length_error("number of samples requested from linspace must be nonzero");
    marray<double,1> linpoints{div};
    double step_lin = (Emax - Emin)/double(div-1);
    
    double EE = Emin;
    for(unsigned int i=0; i<div-1; i++, EE+=step_lin)
        linpoints[i] = EE;
    linpoints[div-1] = Emax;
	
    return linpoints;
}

marray<double,1> logspace(double Emin,double Emax,unsigned int div){
    if(div==0)
        throw std::length_error("number of samples requested from logspace must be nonzero");
    marray<double,1> logpoints{div};
    double Emin_log,Emax_log;
    Emin_log = log(Emin);
    Emax_log = log(Emax);
    
    double step_log = (Emax_log - Emin_log)/double(div-1);
    
    logpoints[0]=Emin;
    double EE = Emin_log+step_log;
    for(unsigned int i=1; i<div-1; i++, EE+=step_log)
        logpoints[i] = exp(EE);
    logpoints[div-1]=Emax;
    return logpoints;
}

// additional GSL-like tools

void gsl_matrix_complex_conjugate(gsl_matrix_complex *cm)
{
  gsl_complex z;
  size_t i, j;
  for (i = 0; i < cm->size1; i++) {
    for (j = 0; j < cm->size2; j++) {
      z = gsl_matrix_complex_get(cm, i, j);
      gsl_matrix_complex_set(cm, i, j, gsl_complex_conjugate(z));
    }
  }
}

void gsl_matrix_complex_print(gsl_matrix_complex* matrix){
    for(unsigned int i = 0; i < matrix->size1; i++){
        for(unsigned int j = 0; j < matrix->size2; j++){
            std::cout << gsl_matrix_complex_get(matrix,i,j).dat[0] <<
            "+i" << gsl_matrix_complex_get(matrix,i,j).dat[1] << " ";
        }
        std::cout << std::endl;
    }
}

void gsl_matrix_complex_change_basis_UMUC(gsl_matrix_complex* U, gsl_matrix_complex* M){
    unsigned int numneu = U->size1;
    gsl_matrix_complex *U1 = gsl_matrix_complex_alloc(numneu,numneu);
    gsl_matrix_complex *U2 = gsl_matrix_complex_alloc(numneu,numneu);
    gsl_matrix_complex_memcpy(U1,U);
    gsl_matrix_complex_memcpy(U2,U);

    gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(numneu,numneu);

    // doing : U M U^dagger

    gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,
                   gsl_complex_rect(1.0,0.0),M,
                   U1,gsl_complex_rect(0.0,0.0),T1);
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,
                   gsl_complex_rect(1.0,0.0),U2,
                   T1,gsl_complex_rect(0.0,0.0),M);
    // now H_current is in the interaction basis of the mass basis

    gsl_matrix_complex_free(U1);
    gsl_matrix_complex_free(U2);
    gsl_matrix_complex_free(T1);
}

void gsl_matrix_complex_change_basis_UCMU(gsl_matrix_complex* U, gsl_matrix_complex* M){
    unsigned int numneu = U->size1;
    gsl_matrix_complex *U1 = gsl_matrix_complex_alloc(numneu,numneu);
    gsl_matrix_complex *U2 = gsl_matrix_complex_alloc(numneu,numneu);
    gsl_matrix_complex_memcpy(U1,U);
    gsl_matrix_complex_memcpy(U2,U);

    gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(numneu,numneu);

    // doing : U M U^dagger

    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,
                   gsl_complex_rect(1.0,0.0),M,
                   U1,gsl_complex_rect(0.0,0.0),T1);
    gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,
                   gsl_complex_rect(1.0,0.0),U2,
                   T1,gsl_complex_rect(0.0,0.0),M);
    // now H_current is in the interaction basis of the mass basis

    gsl_matrix_complex_free(U1);
    gsl_matrix_complex_free(U2);
    gsl_matrix_complex_free(T1);
}
    
//Numerik-Algorithmen: Verfahren, Beispiele, Anwendungen
//Engeln-M{\"u}llges, Gisela and Niederdrenk, Klaus and Wodicka, Reinhard
//https://doi.org/10.1007/978-3-642-13473-9_11
void AkimaSpline::initialize(const double* x, const double* y, const unsigned int n){
    if(n<3)
    	throw std::runtime_error("At least 3 points are required to construct the Akima spline interpolation");
    segments.reserve(n-1);
    segments.resize(0);
    
    //a lambda, basically, but needs to be able to call itself recursively
    struct mHelper{
        const double* x;
        const double* y;
        const unsigned int n;
        double operator()(int i) const{
            const mHelper& m=*this;
            if(i==-2)
                return(3*m(0)-2*m(1));
            if(i==-1)
                return(2*m(0)-m(1));
            if(i==n-1)
                return(2*m(n-2)-m(n-3));
            if(i==n)
                return(3*m(n-2)-2*m(n-3));
            return((y[i+1]-y[i])/(x[i+1]-x[i]));
        }
    } m{x,y,n};
    
    double L_ip1=-1, NE_ip1=-1;
    double m_im2, m_im1, m_i, m_ip1, m_ip2;
    for(int i=0; i<n-1; i++){
        double h_i=x[i+1]-x[i];
        if(h_i<0)
            throw std::runtime_error("AkimaSpline: absissa values not ordered");
        if(h_i==0)
            throw std::runtime_error("AkimaSpline: absissa values not distinct");
        
        //We mostly reuse m values from previous iterations, but the first
        //iteration must start the process, and all subsequent iterations
        //must compute the new m_{i+2}:
        if(i==0){
            m_im2=m(i-2);
            m_im1=m(i-1);
            m_i=m(i);
            m_ip1=m(i+1);
        }
        m_ip2=m(i+2);
        
        double tRi;
        {
            double L_i, NE_i;
            if(L_ip1>=0){ //take values recorded from last iteration
                L_i=L_ip1;
                NE_i=NE_ip1;
            }
            else{ //no previous values available
                L_i=std::abs(m_im2-m_im1);
                NE_i=L_i+std::abs(m_i-m_ip1);
            }
            if(NE_i>0){
                double alpha=L_i/NE_i;
                tRi=(1-alpha)*m_im1+alpha*m_i;
            }
            else
                tRi=m_i;
        }
        double tLip1;
        {
            L_ip1=std::abs(m_im1-m_i);
            NE_ip1=L_ip1+std::abs(m_ip1-m_ip2);
            if(NE_ip1>0){
                double alpha=L_ip1/NE_ip1;
                tLip1=(1-alpha)*m_i+alpha*m_ip1;
            }
            else
                tLip1=m(i);
        }
        
        segment seg;
        seg.x=x[i];
        seg.a0=y[i];
        seg.a1=tRi;
        seg.a2=(3*m_i-2*tRi-tLip1)/h_i;
        seg.a3=(tRi+tLip1-2*m_i)/(h_i*h_i);
        segments.push_back(seg);
        
        //shift values for next iteration
        m_im2=m_im1;
        m_im1=m_i;
        m_i=m_ip1;
        m_ip1=m_ip2;
    }
    xLast=x[n-1];
    yLast=y[n-1];
}

marray<double,1> AkimaSpline::getAbscissas() const{
	marray<double,1> x({segments.size()+1});
	auto it=x.begin();
	for(const auto& s : segments)
		*(it++)=s.x;
	x[segments.size()]=xLast;
	return x;
}

marray<double,1> AkimaSpline::getOrdinates() const{
	marray<double,1> y({segments.size()+1});
	auto it=y.begin();
	for(const auto& s : segments)
		*(it++)=s.a0;
	y[segments.size()]=yLast;
	return y;
}

marray<double,1> AkimaSpline::getAbscissaDerivatives() const{
	marray<double,1> d({segments.size()+1});
	auto it=d.begin();
	for(const auto& s : segments)
		*(it++)=s.a1;
	//this is the only one we don't have precalculated:
	d[segments.size()]=segments.back().derivative(xLast);
	return d;
}

void AkimaSpline::getAbscissaDerivatives(marray<double,1>& d) const{
	d.resize(0,segments.size()+1);
	auto it=d.begin();
	for(const auto& s : segments)
		*(it++)=s.a1;
	//this is the only one we don't have precalculated:
	d[segments.size()]=segments.back().derivative(xLast);
}

BiCubicInterpolator::BiCubicInterpolator(marray<double,2> data_, 
                                         marray<double,1> xcoords_, 
                                         marray<double,1> ycoords_):
data(std::move(data_)),
xcoords(std::move(xcoords_)),
ycoords(std::move(ycoords_)),
dzdx(data.get_extents(),detail::no_init),
dzdy(data.get_extents(),detail::no_init),
dzdxdy(data.get_extents(),detail::no_init)
{	
	if(xcoords.extent(0) != data.extent(1))
		throw std::runtime_error("Number of x-coordinates does not match interpolation grid size");
	if(ycoords.extent(0) != data.extent(0))
		throw std::runtime_error("Number of y-coordinates does not match interpolation grid size");
	
	//compute derivatives at all tabulated points by spline interpolation
	std::size_t iMax=data.extent(1);
	std::size_t jMax=data.extent(0);
	AkimaSpline spline;
	
	marray<double,1> xslice(xcoords.get_extents());
	marray<double,1> yslice(ycoords.get_extents());
	
	for(std::size_t i=0; i!=iMax; i++){
		for(std::size_t j=0; j!=jMax; j++)
			yslice[j]=data[j][i];
		spline.initialize(ycoords.get_data(),yslice.get_data(),jMax);
		spline.getAbscissaDerivatives(yslice);
		for(std::size_t j=0; j!=jMax; j++)
			dzdy[j][i]=yslice[j];
	}
	
	for(std::size_t j=0; j!=jMax; j++){
		for(std::size_t i=0; i!=iMax; i++)
			xslice[i]=data[j][i];
		spline.initialize(xcoords.get_data(),xslice.get_data(),iMax);
		spline.getAbscissaDerivatives(xslice);
		for(std::size_t i=0; i!=iMax; i++)
			dzdx[j][i]=xslice[i];
	}
	
	for(std::size_t j=0; j!=jMax; j++){
		for(std::size_t i=0; i!=iMax; i++)
			xslice[i]=dzdy[j][i];
		spline.initialize(xcoords.get_data(),xslice.get_data(),iMax);
		spline.getAbscissaDerivatives(xslice);
		for(std::size_t i=0; i!=iMax; i++)
			dzdxdy[j][i]=xslice[i];
	}
}

double BiCubicInterpolator::operator()(double x, double y) const{
	if(data.empty())
		return 0;
	
	//figure out which grid cell we are in
	auto x_it=std::upper_bound(xcoords.begin(),xcoords.end(),x);
	std::size_t i=x<xcoords[0] ? 0 : std::distance(xcoords.begin(),x_it)-1;
	if(i==xcoords.extent(0)-1)
		i--;
	auto y_it=std::upper_bound(ycoords.begin(),ycoords.end(),y);
	std::size_t j=y<ycoords[0] ? 0 : std::distance(ycoords.begin(),y_it)-1;
	if(j==ycoords.extent(0)-1)
		j--;
	
	//cell dimensions
	double dx=xcoords[i+1]-xcoords[i];
	double dy=ycoords[j+1]-ycoords[j];
	//compute local coordinates
	double t=(x-xcoords[i])/dx;
	double u=(y-ycoords[j])/dy;
	//Numbering scheme, following NR:
	//   |        |   
	//---3--------2--- y[j+1]
	//   |        |   
	//   |        | dy
	//   |        |   
	//---0--------1--- y[j+2]
	//   |   dx   |   
	//  x[i]    x[i+1]
	double z[4]={
		data[j  ][i  ],
		data[j  ][i+1],
		data[j+1][i+1],
		data[j+1][i  ],
	};
	double dzdt[4]={
		dzdx[j  ][i  ]*dx, // dt/dx=1/dx => dx/dt=dx, dz/dt=dz/dx*dx/dt=dz/dx*dx
		dzdx[j  ][i+1]*dx,
		dzdx[j+1][i+1]*dx,
		dzdx[j+1][i  ]*dx,
	};
	double dzdu[4]={
		dzdy[j  ][i  ]*dy, // du/dy=1/dy => dy/du=dx, dz/du=dz/dy*dy/du=dz/dy*dy
		dzdy[j  ][i+1]*dy,
		dzdy[j+1][i+1]*dy,
		dzdy[j+1][i  ]*dy,
	};
	double dzdtdu[4]={
		dzdxdy[j  ][i  ]*dx*dy, // d2z/dtdu = d2z/dxdy * dx/dt * dy/du = d2z/dxdy*dx*dy
		dzdxdy[j  ][i+1]*dx*dy,
		dzdxdy[j+1][i+1]*dx*dy,
		dzdxdy[j+1][i  ]*dx*dy,
	};
	
	double c0[4]={
		z[0],
		dzdu[0],
		-3*z[0] + 3*z[3] - 2*dzdu[0] - dzdu[3],
		2*z[0] - 2*z[3] + dzdu[0] + dzdu[3],
	};
	double c1[4]={
		dzdt[0],
		dzdtdu[0],
		-3*dzdt[0] + 3*dzdt[3] - 2*dzdtdu[0] - dzdtdu[3],
		2*dzdt[0] - 2*dzdt[3] + dzdtdu[0] + dzdtdu[3],
	};
	double c2[4]={
		-3*z[0] + 3*z[1] - 2*dzdt[0] - dzdt[1],
		-3*dzdu[0] + 3*dzdu[1] - 2*dzdtdu[0] - dzdtdu[1],
		9*z[0] - 9*z[1] + 9*z[2] - 9*z[3] + 6*dzdt[0] + 3*dzdt[1] - 3*dzdt[2] - 6*dzdt[3] + 6*dzdu[0] - 6*dzdu[1] - 3*dzdu[2] + 3*dzdu[3] + 4*dzdtdu[0] + 2*dzdtdu[1] + dzdtdu[2] + 2*dzdtdu[3],
		-6*z[0] + 6*z[1] - 6*z[2] + 6*z[3] - 4*dzdt[0] - 2*dzdt[1] + 2*dzdt[2] + 4*dzdt[3] - 3*dzdu[0] + 3*dzdu[1] + 3*dzdu[2] - 3*dzdu[3] - 2*dzdtdu[0] - dzdtdu[1] - dzdtdu[2] - 2*dzdtdu[3],
	};
	double c3[4]={
		2*z[0] - 2*z[1] + dzdt[0] + dzdt[1],
		2*dzdu[0] - 2*dzdu[1] + dzdtdu[0] + dzdtdu[1],
		-6*z[0] + 6*z[1] - 6*z[2] + 6*z[3] - 3*dzdt[0] - 3*dzdt[1] + 3*dzdt[2] + 3*dzdt[3] - 4*dzdu[0] + 4*dzdu[1] + 2*dzdu[2] - 2*dzdu[3] - 2*dzdtdu[0] - 2*dzdtdu[1] - dzdtdu[2] - dzdtdu[3],
		4*z[0] - 4*z[1] + 4*z[2] - 4*z[3] + 2*dzdt[0] + 2*dzdt[1] - 2*dzdt[2] - 2*dzdt[3] + 2*dzdu[0] - 2*dzdu[1] - 2*dzdu[2] + 2*dzdu[3] + dzdtdu[0] + dzdtdu[1] + dzdtdu[2] + dzdtdu[3],
	};
	
	double result   = ((c3[3]*u+c3[2])*u+c3[1])*u+c3[0];
	result=result*t + ((c2[3]*u+c2[2])*u+c2[1])*u+c2[0];
	result=result*t + ((c1[3]*u+c1[2])*u+c1[1])*u+c1[0];
	result=result*t + ((c0[3]*u+c0[2])*u+c0[1])*u+c0[0];
	return result;
}

TriCubicInterpolator::TriCubicInterpolator(marray<double,3> data_,
                                           marray<double,1> xcoords_,
                                           marray<double,1> ycoords_,
                                           marray<double,1> zcoords_):
xcoords(std::move(xcoords_)),
ycoords(std::move(ycoords_)),
zcoords(std::move(zcoords_)),
data(std::move(data_)),
dfdx(data.get_extents(),nusquids::detail::no_init),
dfdy(data.get_extents(),nusquids::detail::no_init),
dfdz(data.get_extents(),nusquids::detail::no_init),
d2fdxdy(data.get_extents(),nusquids::detail::no_init),
d2fdxdz(data.get_extents(),nusquids::detail::no_init),
d2fdydz(data.get_extents(),nusquids::detail::no_init),
d3fdxdydz(data.get_extents(),nusquids::detail::no_init)
{
	if(xcoords.extent(0) != data.extent(2))
		throw std::runtime_error("Number of x-coordinates does not match interpolation grid size");
	if(ycoords.extent(0) != data.extent(1))
		throw std::runtime_error("Number of y-coordinates does not match interpolation grid size");
	if(zcoords.extent(0) != data.extent(0))
		throw std::runtime_error("Number of z-coordinates does not match interpolation grid size");
	
	std::size_t iMax=data.extent(2);
	std::size_t jMax=data.extent(1);
	std::size_t kMax=data.extent(0);
	AkimaSpline spline;
	marray<double,1> xslice(xcoords.get_extents());
	marray<double,1> yslice(ycoords.get_extents());
	marray<double,1> zslice(zcoords.get_extents());
	
	for(std::size_t k=0; k!=kMax; k++){
		for(std::size_t j=0; j!=jMax; j++){
			for(std::size_t i=0; i!=iMax; i++)
				xslice[i]=data[k][j][i];
			spline.initialize(xcoords.get_data(),xslice.get_data(),iMax);
			spline.getAbscissaDerivatives(xslice);
			for(std::size_t i=0; i!=iMax; i++)
				dfdx[k][j][i]=xslice[i];
		}
	}
	
	for(std::size_t k=0; k!=kMax; k++){
		for(std::size_t i=0; i!=iMax; i++){
			for(std::size_t j=0; j!=jMax; j++)
				yslice[j]=data[k][j][i];
			spline.initialize(ycoords.get_data(),yslice.get_data(),jMax);
			spline.getAbscissaDerivatives(yslice);
			for(std::size_t j=0; j!=jMax; j++)
				dfdy[k][j][i]=yslice[j];
		}
	}
	
	for(std::size_t j=0; j!=jMax; j++){
		for(std::size_t i=0; i!=iMax; i++){
			for(std::size_t k=0; k!=kMax; k++)
				zslice[k]=data[k][j][i];
			spline.initialize(zcoords.get_data(),zslice.get_data(),kMax);
			spline.getAbscissaDerivatives(zslice);
			for(std::size_t k=0; k!=kMax; k++)
				dfdz[k][j][i]=zslice[k];
		}
	}
	
	for(std::size_t k=0; k!=kMax; k++){
		for(std::size_t j=0; j!=jMax; j++){
			for(std::size_t i=0; i!=iMax; i++)
				xslice[i]=dfdy[k][j][i];
			spline.initialize(xcoords.get_data(),xslice.get_data(),iMax);
			spline.getAbscissaDerivatives(xslice);
			for(std::size_t i=0; i!=iMax; i++)
				d2fdxdy[k][j][i]=xslice[i];
		}
	}
	
	for(std::size_t k=0; k!=kMax; k++){
		for(std::size_t j=0; j!=jMax; j++){
			for(std::size_t i=0; i!=iMax; i++)
				xslice[i]=dfdz[k][j][i];
			spline.initialize(xcoords.get_data(),xslice.get_data(),iMax);
			spline.getAbscissaDerivatives(xslice);
			for(std::size_t i=0; i!=iMax; i++)
				d2fdxdz[k][j][i]=xslice[i];
		}
	}
	
	for(std::size_t k=0; k!=kMax; k++){
		for(std::size_t i=0; i!=iMax; i++){			
			for(std::size_t j=0; j!=jMax; j++)
				yslice[j]=dfdz[k][j][i];
			spline.initialize(ycoords.get_data(),yslice.get_data(),jMax);
			spline.getAbscissaDerivatives(yslice);
			for(std::size_t j=0; j!=jMax; j++)
				d2fdydz[k][j][i]=xslice[j];
		}
	}
	
	for(std::size_t k=0; k!=kMax; k++){
		for(std::size_t j=0; j!=jMax; j++){
			for(std::size_t i=0; i!=iMax; i++)
				xslice[i]=d2fdydz[k][j][i];
			spline.initialize(xcoords.get_data(),xslice.get_data(),iMax);
			spline.getAbscissaDerivatives(xslice);
			for(std::size_t i=0; i!=iMax; i++)
				d3fdxdydz[k][j][i]=xslice[i];
		}
	}
}

void TriCubicInterpolator::collect_data(double d[64], std::size_t i, std::size_t j, std::size_t k, double dx, double dy, double dz) const{
	//libtricubic vertex numbering scheme, which is *not* an extension of NR's 2D scheme
	//    z ^
	//      |
	//     4|          6
	//      o---------o
	//     /|        /|
	//  5 / |     7 / | 
	//   o---------o  |
	//   |  |      |  |
	//   |  o------|--o-->
	//   | / 0     | / 2  y
	//   |/        |/
	//   o---------o
	//  / 1         3
	// L x
	
	d[ 0]=data[k  ][j  ][i  ];
	d[ 1]=data[k  ][j  ][i+1];
	d[ 2]=data[k  ][j+1][i  ];
	d[ 3]=data[k  ][j+1][i+1];
	d[ 4]=data[k+1][j  ][i  ];
	d[ 5]=data[k+1][j  ][i+1];
	d[ 6]=data[k+1][j+1][i  ];
	d[ 7]=data[k+1][j+1][i+1];
	d[ 8]=dfdx[k  ][j  ][i  ]*dx;
	d[ 9]=dfdx[k  ][j  ][i+1]*dx;
	d[10]=dfdx[k  ][j+1][i  ]*dx;
	d[11]=dfdx[k  ][j+1][i+1]*dx;
	d[12]=dfdx[k+1][j  ][i  ]*dx;
	d[13]=dfdx[k+1][j  ][i+1]*dx;
	d[14]=dfdx[k+1][j+1][i  ]*dx;
	d[15]=dfdx[k+1][j+1][i+1]*dx;
	d[16]=dfdy[k  ][j  ][i  ]*dy;
	d[17]=dfdy[k  ][j  ][i+1]*dy;
	d[18]=dfdy[k  ][j+1][i  ]*dy;
	d[19]=dfdy[k  ][j+1][i+1]*dy;
	d[20]=dfdy[k+1][j  ][i  ]*dy;
	d[21]=dfdy[k+1][j  ][i+1]*dy;
	d[22]=dfdy[k+1][j+1][i  ]*dy;
	d[23]=dfdy[k+1][j+1][i+1]*dy;
	d[24]=dfdz[k  ][j  ][i  ]*dz;
	d[25]=dfdz[k  ][j  ][i+1]*dz;
	d[26]=dfdz[k  ][j+1][i  ]*dz;
	d[27]=dfdz[k  ][j+1][i+1]*dz;
	d[28]=dfdz[k+1][j  ][i  ]*dz;
	d[29]=dfdz[k+1][j  ][i+1]*dz;
	d[30]=dfdz[k+1][j+1][i  ]*dz;
	d[31]=dfdz[k+1][j+1][i+1]*dz;
	d[32]=d2fdxdy[k  ][j  ][i  ]*dx*dy;
	d[33]=d2fdxdy[k  ][j  ][i+1]*dx*dy;
	d[34]=d2fdxdy[k  ][j+1][i  ]*dx*dy;
	d[35]=d2fdxdy[k  ][j+1][i+1]*dx*dy;
	d[36]=d2fdxdy[k+1][j  ][i  ]*dx*dy;
	d[37]=d2fdxdy[k+1][j  ][i+1]*dx*dy;
	d[38]=d2fdxdy[k+1][j+1][i  ]*dx*dy;
	d[39]=d2fdxdy[k+1][j+1][i+1]*dx*dy;
	d[40]=d2fdxdz[k  ][j  ][i  ]*dx*dz;
	d[41]=d2fdxdz[k  ][j  ][i+1]*dx*dz;
	d[42]=d2fdxdz[k  ][j+1][i  ]*dx*dz;
	d[43]=d2fdxdz[k  ][j+1][i+1]*dx*dz;
	d[44]=d2fdxdz[k+1][j  ][i  ]*dx*dz;
	d[45]=d2fdxdz[k+1][j  ][i+1]*dx*dz;
	d[46]=d2fdxdz[k+1][j+1][i  ]*dx*dz;
	d[47]=d2fdxdz[k+1][j+1][i+1]*dx*dz;
	d[48]=d2fdydz[k  ][j  ][i  ]*dy*dz;
	d[49]=d2fdydz[k  ][j  ][i+1]*dy*dz;
	d[50]=d2fdydz[k  ][j+1][i  ]*dy*dz;
	d[51]=d2fdydz[k  ][j+1][i+1]*dy*dz;
	d[52]=d2fdydz[k+1][j  ][i  ]*dy*dz;
	d[53]=d2fdydz[k+1][j  ][i+1]*dy*dz;
	d[54]=d2fdydz[k+1][j+1][i  ]*dy*dz;
	d[55]=d2fdydz[k+1][j+1][i+1]*dy*dz;
	d[56]=d3fdxdydz[k  ][j  ][i  ]*dx*dy*dz;
	d[57]=d3fdxdydz[k  ][j  ][i+1]*dx*dy*dz;
	d[58]=d3fdxdydz[k  ][j+1][i  ]*dx*dy*dz;
	d[59]=d3fdxdydz[k  ][j+1][i+1]*dx*dy*dz;
	d[60]=d3fdxdydz[k+1][j  ][i  ]*dx*dy*dz;
	d[61]=d3fdxdydz[k+1][j  ][i+1]*dx*dy*dz;
	d[62]=d3fdxdydz[k+1][j+1][i  ]*dx*dy*dz;
	d[63]=d3fdxdydz[k+1][j+1][i+1]*dx*dy*dz;
}

void TriCubicInterpolator::get_coeffs(double c[64], double d[64]){
	c[ 0] = d[0];
	c[ 1] = d[8];
	c[ 2] = 3*(-d[0] + d[1]) - 2*d[8] - d[9];
	c[ 3] = 2*(d[0] - d[1]) + d[8] + d[9];
	c[ 4] = d[16];
	c[ 5] = d[32];
	c[ 6] = 3*(-d[16] + d[17]) - 2*d[32] - d[33];
	c[ 7] = 2*(d[16] - d[17]) + d[32] + d[33];
	c[ 8] = 3*(-d[0] + d[2]) - 2*d[16] - d[18];
	c[ 9] = 3*(-d[8] + d[10]) - 2*d[32] - d[34];
	c[10] = 9*(d[0] - d[1] - d[2] + d[3]) + 6*(d[8] - d[10] + d[16] - d[17]) + 4*d[32] + 3*(d[9] - d[11] + d[18] - d[19]) + 2*(d[33] + d[34]) + d[35];
	c[11] = 6*(-d[0] + d[1] + d[2] - d[3]) + 4*(-d[16] + d[17]) + 3*(-d[8] - d[9] + d[10] + d[11]) + 2*(-d[18] + d[19] - d[32] - d[33]) - d[34] - d[35];
	c[12] = 2*(d[0] - d[2]) + d[16] + d[18];
	c[13] = 2*(d[8] - d[10]) + d[32] + d[34];
	c[14] = 6*(-d[0] + d[1] + d[2] - d[3]) + 4*(-d[8] + d[10]) + 3*(-d[16] + d[17] - d[18] + d[19]) + 2*(-d[9] + d[11] - d[32] - d[34]) - d[33] - d[35];
	c[15] = 4*(d[0] - d[1] - d[2] + d[3]) + 2*(d[8] + d[9] - d[10] - d[11] + d[16] - d[17] + d[18] - d[19]) + d[32] + d[33] + d[34] + d[35];
	c[16] = d[24];
	c[17] = d[40];
	c[18] = 3*(-d[24] + d[25]) - 2*d[40] - d[41];
	c[19] = 2*(d[24] - d[25]) + d[40] + d[41];
	c[20] = d[48];
	c[21] = d[56];
	c[22] = 3*(-d[48] + d[49]) - 2*d[56] - d[57];
	c[23] = 2*(d[48] - d[49]) + d[56] + d[57];
	c[24] = 3*(-d[24] + d[26]) - 2*d[48] - d[50];
	c[25] = 3*(-d[40] + d[42]) - 2*d[56] - d[58];
	c[26] = 9*(d[24] - d[25] - d[26] + d[27]) + 6*(d[40] - d[42] + d[48] - d[49]) + 4*d[56] + 3*(d[41] - d[43] + d[50] - d[51]) + 2*(d[57] + d[58]) + d[59];
	c[27] = 6*(-d[24] + d[25] + d[26] - d[27]) + 4*(-d[48] + d[49]) + 3*(-d[40] - d[41] + d[42] + d[43]) + 2*(-d[50] + d[51] - d[56] - d[57]) - d[58] - d[59];
	c[28] = 2*(d[24] - d[26]) + d[48] + d[50];
	c[29] = 2*(d[40] - d[42]) + d[56] + d[58];
	c[30] = 6*(-d[24] + d[25] + d[26] - d[27]) + 4*(-d[40] + d[42]) + 3*(-d[48] + d[49] - d[50] + d[51]) + 2*(-d[41] + d[43] - d[56] - d[58]) - d[57] - d[59];
	c[31] = 4*(d[24] - d[25] - d[26] + d[27]) + 2*(d[40] + d[41] - d[42] - d[43] + d[48] - d[49] + d[50] - d[51]) + d[56] + d[57] + d[58] + d[59];
	c[32] = 3*(-d[0] + d[4]) - 2*d[24] - d[28];
	c[33] = 3*(-d[8] + d[12]) - 2*d[40] - d[44];
	c[34] = 9*(d[0] - d[1] - d[4] + d[5]) + 6*(d[8] - d[12] + d[24] - d[25]) + 4*d[40] + 3*(d[9] - d[13] + d[28] - d[29]) + 2*(d[41] + d[44]) + d[45];
	c[35] = 6*(-d[0] + d[1] + d[4] - d[5]) + 4*(-d[24] + d[25]) + 3*(-d[8] - d[9] + d[12] + d[13]) + 2*(-d[28] + d[29] - d[40] - d[41]) - d[44] - d[45];
	c[36] = 3*(-d[16] + d[20]) - 2*d[48] - d[52];
	c[37] = 3*(-d[32] + d[36]) - 2*d[56] - d[60];
	c[38] = 9*(d[16] - d[17] - d[20] + d[21]) + 6*(d[32] - d[36] + d[48] - d[49]) + 4*d[56] + 3*(d[33] - d[37] + d[52] - d[53]) + 2*(d[57] + d[60]) + d[61];
	c[39] = 6*(-d[16] + d[17] + d[20] - d[21]) + 4*(-d[48] + d[49]) + 3*(-d[32] - d[33] + d[36] + d[37]) + 2*(-d[52] + d[53] - d[56] - d[57]) - d[60] - d[61];
	c[40] = 9*(d[0] - d[2] - d[4] + d[6]) + 6*(d[16] - d[20] + d[24] - d[26]) + 4*d[48] + 3*(d[18] - d[22] + d[28] - d[30]) + 2*(d[50] + d[52]) + d[54];
	c[41] = 9*(d[8] - d[10] - d[12] + d[14]) + 6*(d[32] - d[36] + d[40] - d[42]) + 4*d[56] + 3*(d[34] - d[38] + d[44] - d[46]) + 2*(d[58] + d[60]) + d[62];
	c[42] = 27*(-d[0] + d[1] + d[2] - d[3] + d[4] - d[5] - d[6] + d[7]) + 18*(-d[8] + d[10] + d[12] - d[14] - d[16] + d[17] + d[20] - d[21] - d[24] + d[25] + d[26] - d[27]) + 12*(-d[32] + d[36] - d[40] + d[42] - d[48] + d[49]) + 9*(-d[9] + d[11] + d[13] - d[15] - d[18] + d[19] + d[22] - d[23] - d[28] + d[29] + d[30] - d[31]) - 8*d[56] + 6*(-d[33] - d[34] + d[37] + d[38] - d[41] + d[43] - d[44] + d[46] - d[50] + d[51] - d[52] + d[53]) + 4*(-d[57] - d[58] - d[60]) + 3*(-d[35] + d[39] - d[45] + d[47] - d[54] + d[55]) + 2*(-d[59] - d[61] - d[62]) - d[63];
	c[43] = 18*(d[0] - d[1] - d[2] + d[3] - d[4] + d[5] + d[6] - d[7]) + 12*(d[16] - d[17] - d[20] + d[21] + d[24] - d[25] - d[26] + d[27]) + 9*(d[8] + d[9] - d[10] - d[11] - d[12] - d[13] + d[14] + d[15]) + 8*(d[48] - d[49]) + 6*(d[18] - d[19] - d[22] + d[23] + d[28] - d[29] - d[30] + d[31] + d[32] + d[33] - d[36] - d[37] + d[40] + d[41] - d[42] - d[43]) + 4*(d[50] - d[51] + d[52] - d[53] + d[56] + d[57]) + 3*(d[34] + d[35] - d[38] - d[39] + d[44] + d[45] - d[46] - d[47]) + 2*(d[54] - d[55] + d[58] + d[59] + d[60] + d[61]) + d[62] + d[63];
	c[44] = 6*(-d[0] + d[2] + d[4] - d[6]) + 4*(-d[24] + d[26]) + 3*(-d[16] - d[18] + d[20] + d[22]) + 2*(-d[28] + d[30] - d[48] - d[50]) - d[52] - d[54];
	c[45] = 6*(-d[8] + d[10] + d[12] - d[14]) + 4*(-d[40] + d[42]) + 3*(-d[32] - d[34] + d[36] + d[38]) + 2*(-d[44] + d[46] - d[56] - d[58]) - d[60] - d[62];
	c[46] = 18*(d[0] - d[1] - d[2] + d[3] - d[4] + d[5] + d[6] - d[7]) + 12*(d[8] - d[10] - d[12] + d[14] + d[24] - d[25] - d[26] + d[27]) + 9*(d[16] - d[17] + d[18] - d[19] - d[20] + d[21] - d[22] + d[23]) + 8*(d[40] - d[42]) + 6*(d[9] - d[11] - d[13] + d[15] + d[28] - d[29] - d[30] + d[31] + d[32] + d[34] - d[36] - d[38] + d[48] - d[49] + d[50] - d[51]) + 4*(d[41] - d[43] + d[44] - d[46] + d[56] + d[58]) + 3*(d[33] + d[35] - d[37] - d[39] + d[52] - d[53] + d[54] - d[55]) + 2*(d[45] - d[47] + d[57] + d[59] + d[60] + d[62]) + d[61] + d[63];
	c[47] = 12*(-d[0] + d[1] + d[2] - d[3] + d[4] - d[5] - d[6] + d[7]) + 8*(-d[24] + d[25] + d[26] - d[27]) + 6*(-d[8] - d[9] + d[10] + d[11] + d[12] + d[13] - d[14] - d[15] - d[16] + d[17] - d[18] + d[19] + d[20] - d[21] + d[22] - d[23]) + 4*(-d[28] + d[29] + d[30] - d[31] - d[40] - d[41] + d[42] + d[43] - d[48] + d[49] - d[50] + d[51]) + 3*(-d[32] - d[33] - d[34] - d[35] + d[36] + d[37] + d[38] + d[39]) + 2*(-d[44] - d[45] + d[46] + d[47] - d[52] + d[53] - d[54] + d[55] - d[56] - d[57] - d[58] - d[59]) - d[60] - d[61] - d[62] - d[63];
	c[48] = 2*(d[0] - d[4]) + d[24] + d[28];
	c[49] = 2*(d[8] - d[12]) + d[40] + d[44];
	c[50] = 6*(-d[0] + d[1] + d[4] - d[5]) + 4*(-d[8] + d[12]) + 3*(-d[24] + d[25] - d[28] + d[29]) + 2*(-d[9] + d[13] - d[40] - d[44]) - d[41] - d[45];
	c[51] = 4*(d[0] - d[1] - d[4] + d[5]) + 2*(d[8] + d[9] - d[12] - d[13] + d[24] - d[25] + d[28] - d[29]) + d[40] + d[41] + d[44] + d[45];
	c[52] = 2*(d[16] - d[20]) + d[48] + d[52];
	c[53] = 2*(d[32] - d[36]) + d[56] + d[60];
	c[54] = 6*(-d[16] + d[17] + d[20] - d[21]) + 4*(-d[32] + d[36]) + 3*(-d[48] + d[49] - d[52] + d[53]) + 2*(-d[33] + d[37] - d[56] - d[60]) - d[57] - d[61];
	c[55] = 4*(d[16] - d[17] - d[20] + d[21]) + 2*(d[32] + d[33] - d[36] - d[37] + d[48] - d[49] + d[52] - d[53]) + d[56] + d[57] + d[60] + d[61];
	c[56] = 6*(-d[0] + d[2] + d[4] - d[6]) + 4*(-d[16] + d[20]) + 3*(-d[24] + d[26] - d[28] + d[30]) + 2*(-d[18] + d[22] - d[48] - d[52]) - d[50] - d[54];
	c[57] = 6*(-d[8] + d[10] + d[12] - d[14]) + 4*(-d[32] + d[36]) + 3*(-d[40] + d[42] - d[44] + d[46]) + 2*(-d[34] + d[38] - d[56] - d[60]) - d[58] - d[62];
	c[58] = 18*(d[0] - d[1] - d[2] + d[3] - d[4] + d[5] + d[6] - d[7]) + 12*(d[8] - d[10] - d[12] + d[14] + d[16] - d[17] - d[20] + d[21]) + 9*(d[24] - d[25] - d[26] + d[27] + d[28] - d[29] - d[30] + d[31]) + 8*(d[32] - d[36]) + 6*(d[9] - d[11] - d[13] + d[15] + d[18] - d[19] - d[22] + d[23] + d[40] - d[42] + d[44] - d[46] + d[48] - d[49] + d[52] - d[53]) + 4*(d[33] + d[34] - d[37] - d[38] + d[56] + d[60]) + 3*(d[41] - d[43] + d[45] - d[47] + d[50] - d[51] + d[54] - d[55]) + 2*(d[35] - d[39] + d[57] + d[58] + d[61] + d[62]) + d[59] + d[63];
	c[59] = 12*(-d[0] + d[1] + d[2] - d[3] + d[4] - d[5] - d[6] + d[7]) + 8*(-d[16] + d[17] + d[20] - d[21]) + 6*(-d[8] - d[9] + d[10] + d[11] + d[12] + d[13] - d[14] - d[15] - d[24] + d[25] + d[26] - d[27] - d[28] + d[29] + d[30] - d[31]) + 4*(-d[18] + d[19] + d[22] - d[23] - d[32] - d[33] + d[36] + d[37] - d[48] + d[49] - d[52] + d[53]) + 3*(-d[40] - d[41] + d[42] + d[43] - d[44] - d[45] + d[46] + d[47]) + 2*(-d[34] - d[35] + d[38] + d[39] - d[50] + d[51] - d[54] + d[55] - d[56] - d[57] - d[60] - d[61]) - d[58] - d[59] - d[62] - d[63];
	c[60] = 4*(d[0] - d[2] - d[4] + d[6]) + 2*(d[16] + d[18] - d[20] - d[22] + d[24] - d[26] + d[28] - d[30]) + d[48] + d[50] + d[52] + d[54];
	c[61] = 4*(d[8] - d[10] - d[12] + d[14]) + 2*(d[32] + d[34] - d[36] - d[38] + d[40] - d[42] + d[44] - d[46]) + d[56] + d[58] + d[60] + d[62];
	c[62] = 12*(-d[0] + d[1] + d[2] - d[3] + d[4] - d[5] - d[6] + d[7]) + 8*(-d[8] + d[10] + d[12] - d[14]) + 6*(-d[16] + d[17] - d[18] + d[19] + d[20] - d[21] + d[22] - d[23] - d[24] + d[25] + d[26] - d[27] - d[28] + d[29] + d[30] - d[31]) + 4*(-d[9] + d[11] + d[13] - d[15] - d[32] - d[34] + d[36] + d[38] - d[40] + d[42] - d[44] + d[46]) + 3*(-d[48] + d[49] - d[50] + d[51] - d[52] + d[53] - d[54] + d[55]) + 2*(-d[33] - d[35] + d[37] + d[39] - d[41] + d[43] - d[45] + d[47] - d[56] - d[58] - d[60] - d[62]) - d[57] - d[59] - d[61] - d[63];
	c[63] = 8*(d[0] - d[1] - d[2] + d[3] - d[4] + d[5] + d[6] - d[7]) + 4*(d[8] + d[9] - d[10] - d[11] - d[12] - d[13] + d[14] + d[15] + d[16] - d[17] + d[18] - d[19] - d[20] + d[21] - d[22] + d[23] + d[24] - d[25] - d[26] + d[27] + d[28] - d[29] - d[30] + d[31]) + 2*(d[32] + d[33] + d[34] + d[35] - d[36] - d[37] - d[38] - d[39] + d[40] + d[41] - d[42] - d[43] + d[44] + d[45] - d[46] - d[47] + d[48] - d[49] + d[50] - d[51] + d[52] - d[53] + d[54] - d[55]) + d[56] + d[57] + d[58] + d[59] + d[60] + d[61] + d[62] + d[63];
}

double TriCubicInterpolator::evaluate_polynomial(double c[64], double t, double u, double v){
	//Empirically, this fully unrolled and factorized form for the evaulation is faster than
	//the deeply nested loop used by libtricubic, or the more optimized form thereof used by
	//https://github.com/igmhub/likely . The code size is large, but it contains fewer floating-
	//point operations (particularly multiplies), has no integer loop book-keeping, and should
	//give significant opportunities to exploit instruction-level parallelism.
	
	//16 polynomials in the x (t) diection
	double p[16];
	p[0]=((c[3]*t+c[2])*t+c[1])*t+c[0];
	p[1]=((c[7]*t+c[6])*t+c[5])*t+c[4];
	p[2]=((c[11]*t+c[10])*t+c[9])*t+c[8];
	p[3]=((c[15]*t+c[14])*t+c[13])*t+c[12];
	p[4]=((c[19]*t+c[18])*t+c[17])*t+c[16];
	p[5]=((c[23]*t+c[22])*t+c[21])*t+c[20];
	p[6]=((c[27]*t+c[26])*t+c[25])*t+c[24];
	p[7]=((c[31]*t+c[30])*t+c[29])*t+c[28];
	p[8]=((c[35]*t+c[34])*t+c[33])*t+c[32];
	p[9]=((c[39]*t+c[38])*t+c[37])*t+c[36];
	p[10]=((c[43]*t+c[42])*t+c[41])*t+c[40];
	p[11]=((c[47]*t+c[46])*t+c[45])*t+c[44];
	p[12]=((c[51]*t+c[50])*t+c[49])*t+c[48];
	p[13]=((c[55]*t+c[54])*t+c[53])*t+c[52];
	p[14]=((c[59]*t+c[58])*t+c[57])*t+c[56];
	p[15]=((c[63]*t+c[62])*t+c[61])*t+c[60];
	//4 polynomials in the y (u) direction
	double q[4];
	q[0]=((p[3]*u+p[2])*u+p[1])*u+p[0];
	q[1]=((p[7]*u+p[6])*u+p[5])*u+p[4];
	q[2]=((p[11]*u+p[10])*u+p[9])*u+p[8];
	q[3]=((p[15]*u+p[14])*u+p[13])*u+p[12];
	//one polynomial in the z (v) direction
	double r=((q[3]*v+q[2])*v+q[1])*v+q[0];
	return r;
}

double TriCubicInterpolator::operator()(double x, double y, double z) const{
	//figure out which grid cell we are in
	auto x_it=std::upper_bound(xcoords.begin(),xcoords.end(),x);
	std::size_t i=x<xcoords[0] ? 0 : std::distance(xcoords.begin(),x_it)-1;
	if(i==xcoords.extent(0)-1)
		i--;
	auto y_it=std::upper_bound(ycoords.begin(),ycoords.end(),y);
	std::size_t j=y<ycoords[0] ? 0 : std::distance(ycoords.begin(),y_it)-1;
	if(j==ycoords.extent(0)-1)
		j--;
	auto z_it=std::upper_bound(zcoords.begin(),zcoords.end(),z);
	std::size_t k=z<zcoords[0] ? 0 : std::distance(zcoords.begin(),z_it)-1;
	if(k==zcoords.extent(0)-1)
		k--;
	//cell dimensions
	double dx=xcoords[i+1]-xcoords[i];
	double dy=ycoords[j+1]-ycoords[j];
	double dz=zcoords[k+1]-zcoords[k];
	//compute local coordinates
	double t=(x-xcoords[i])/dx;
	double u=(y-ycoords[j])/dy;
	double v=(z-zcoords[k])/dz;
	
	double d[64];
	collect_data(d,i,j,k,dx,dy,dz);
	double c[64];
	get_coeffs(c,d);
	return evaluate_polynomial(c,t,u,v);
}

bool h5ObjectExists(hid_t loc_id, const char* name){
    return H5Lexists(loc_id,name,H5P_DEFAULT)>0 && H5Oexists_by_name(loc_id,name,H5P_DEFAULT)>0;
};
    
template<>
void addH5Attribute<std::string>(hid_t object, std::string name, const std::string& contents){
    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, contents.size());
    hsize_t dim=1;
    hid_t dataspace_id = H5Screate_simple(1, &dim, nullptr);
    hid_t attribute_id = H5Acreate(object,name.c_str(),strtype,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
    H5Awrite(attribute_id, strtype, &contents[0]);
    H5Aclose(attribute_id);
    H5Sclose(dataspace_id);
};

template<>
void readH5Attribute<std::string>(hid_t object, std::string name, std::string& dest){
    hid_t attribute_id = H5Aopen(object, name.c_str(), H5P_DEFAULT);
    hid_t actualType = H5Aget_type(attribute_id);
    H5T_class_t typeClass = H5Tget_class(actualType);
    if(typeClass!=H5T_STRING){
        H5Aclose(attribute_id);
        H5Tclose(actualType);
        throw std::runtime_error("Expected and actual data types for attribute '"+name+"' do not match");
    }
    
    size_t size = H5Tget_size(actualType);
    dest.resize(size);
    
    if(H5Aread(attribute_id, actualType, &dest[0])<0){
        H5Aclose(attribute_id);
        H5Tclose(actualType);
        throw std::runtime_error("Failed to read attribute '"+name+"'");
    }
    
    H5Aclose(attribute_id);
    H5Tclose(actualType);
}

} // close namespace
