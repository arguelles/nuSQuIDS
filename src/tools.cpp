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
