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


#ifndef NUSQUIDS_TOOLS_H
#define NUSQUIDS_TOOLS_H

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <H5Apublic.h>
#include <H5Dpublic.h>
#include <H5Fpublic.h>
#include <H5Ipublic.h>
#include <H5Ppublic.h>
#include <H5Spublic.h>
#include <H5Tpublic.h>
#include "nuSQuIDS/marray.h"

namespace nusquids{

// file read and write
/// \brief Checks if a file exist..
/// @param filename File which exist to check.
bool fexists(const std::string filename);

/// \brief Reads and return the values from a file as a bidimensional array.
/// @param filename Filename to read.
marray<double,2> quickread(const std::string filename);

/// \brief Writes a bidimensional array onto a file.
/// @param filename Filename to write onto.
/// @param array Array to write onto file.
int quickwrite(const std::string filename,const marray<double,2>& array);

/// \brief Construct a linear space
/// @param min Minimum value in the linear span.
/// @param max Maximum value in the linear span.
/// @param samples Number of samples to generate.
marray<double,1> linspace(double min,double max,unsigned int samples);

/// \brief Construct a logarithmic space.
/// @param min Minimum value in the logarithmic span.
/// @param max Maximum value in the logarithmic span.
/// @param samples Number of samples to generate.
marray<double,1> logspace(double min,double max,unsigned int samples);

/// \brief Calculate the complex conjugate of a matrix.
/// @param matrix Matrix to conjugate.
void gsl_matrix_complex_conjugate(gsl_matrix_complex* matrix);

/// \brief Print to std::cout a gsl complex matrix.
/// @param matrix Matrix to print.
void gsl_matrix_complex_print(gsl_matrix_complex* matrix);

/// \brief Perform a unitary rotation.
/// @param U Rotation to apply.
/// @param matrix Matrix to print.
void gsl_matrix_complex_change_basis_UMUC(gsl_matrix_complex* U, gsl_matrix_complex* matrix);

/// \brief Perform a unitary rotation.
/// @param U Rotation to apply.
/// @param matrix Matrix to print.
void gsl_matrix_complex_change_basis_UCMU(gsl_matrix_complex* U, gsl_matrix_complex* matrix);

///\class
///\brief Container for GSL workspace to be use with the integrators.
class IntegrateWorkspace {
private:

public:
    gsl_integration_workspace* ws;
    IntegrateWorkspace(size_t limit) {
        ws=gsl_integration_workspace_alloc(limit);
    }
    ~IntegrateWorkspace() {
        gsl_integration_workspace_free(ws);
    }
};

///\brief One dimensional integral using GSL.
/// @param ws GSL integration workspace.
/// @param f Function to integrate.
/// @param a Lower integration limit.
/// @param b Upper integration limit.
/// @param acc Accuracy parameter.
/// @param max_iter Maximum number of iterations to perform the integral.
template<typename FunctionType>
double integrate(IntegrateWorkspace& ws, FunctionType f, double a, double b, double acc=1e-7, unsigned int max_iter=5000){
    double (*wrapper)(double,void*)=[](double x, void* params){
        FunctionType& f=*static_cast<FunctionType*>(params);
        return(f(x));
    };

    double result, error;
    gsl_function F;
    F.function = wrapper;
    F.params = &f;

    gsl_integration_qags(&F, a, b, 0, acc, max_iter, ws.ws, &result, &error);

    return(result);
}

///\brief One dimensional integral using GSL.
/// @param f Function to integrate.
/// @param a Lower integration limit.
/// @param b Upper integration limit.
/// @param acc Accuracy parameter.
/// @param max_iter Maximum number of iterations to perform the integral.
template<typename FunctionType>
double integrate(FunctionType f, double a, double b, double acc=1e-7, unsigned int max_iter=5000, size_t memory_alloc = 5000){
    IntegrateWorkspace ws(memory_alloc);
    return integrate(ws, f, a, b, acc, max_iter);
}

///\class
///\brief Akima spline implementation
class AkimaSpline{
private:
    struct segment{
        double x;
        double a0, a1, a2, a3;
        
        double operator()(double x) const{
            x-=this->x;
            return((((a3*x)+a2)*x+a1)*x+a0);
        }
        double derivative(double x) const{
            x-=this->x;
            return(((3*a3)*x+2*a2)*x+a1);
        }
        operator double() const{ return(x); }
    };
    std::vector<segment> segments;
    double xLast, yLast;
    template<typename T=void>
    T emitError(const std::string& msg){ throw std::runtime_error(msg); }
public:
	void initialize(const double* x, const double* y, const unsigned int n);
    AkimaSpline(){}
    AkimaSpline(const double* x, const double* y, const unsigned int n){
    	initialize(x,y,n);
    }
    AkimaSpline(const std::vector<double>& x, const std::vector<double>& y):
    AkimaSpline(x.data(),y.data(),
                (x.size()==y.size()?x.size():emitError<std::size_t>("Data array sizes must match"))){}
    
    double operator()(double x) const{
        if(x<segments.front().x)
            return(segments.front()(x));
        auto seg_it=std::upper_bound(segments.begin(),segments.end(),x);
        seg_it--;
        return((*seg_it)(x));
    }
    marray<double,1> getAbscissas() const;
    marray<double,1> getOrdinates() const;
    marray<double,1> getAbscissaDerivatives() const;
    void getAbscissaDerivatives(marray<double,1>& d) const;
};

///\brief A bicubic interpolator for two-dimensional data. 
///One-dimensional Akima splines are used to obtain all necessary derivatives at 
///the grid points. 
struct BiCubicInterpolator{
public:
	BiCubicInterpolator(){}

	BiCubicInterpolator(marray<double,2> data_, 
	                    marray<double,1> xcoords_, 
	                    marray<double,1> ycoords_);
	
	///Compute the interpolated value at the given coordinates. May not return 
	///useful results outside the domain of the interpolation. 
	double operator()(double x, double y) const;
	
	const marray<double,2>& getData() const{ return data; }
	const marray<double,1>& getXCoords() const{ return xcoords; }
	const marray<double,1>& getYCoords() const{ return ycoords; }
private:
	marray<double,2> data;
	marray<double,1> xcoords;
	marray<double,1> ycoords;
	//derivatives
	marray<double,2> dzdx, dzdy, dzdxdy;
};

///\brief A tricubic interpolator for three-dimensional data. 
///One-dimensional Akima splines are used to obtain all necessary derivatives at 
///the grid points. 
class TriCubicInterpolator{
public:
	TriCubicInterpolator();
	///Contruct an interpolator over existing function data for a given grid.
	///\param data Values of the function to be interpolated, evaluated on the grid defined by the
	///            tensor product of the corrdinate arrays. Note that this should be in row-major
	///            order, so that the first index corresponds to the z coordinate, the second to y, 
	///            and the third to x.
	///\param xcoords Values for the x coordinates of the grid points used. THese must be in sorted
	///               order and unique. 
	///\param ycoords Values for the y coordinates of the grid points used. THese must be in sorted
	///               order and unique. 
	///\param zcoords Values for the z coordinates of the grid points used. THese must be in sorted
	///               order and unique. 
	TriCubicInterpolator(marray<double,3> data,
                         marray<double,1> xcoords,
                         marray<double,1> ycoords,
                         marray<double,1> zcoords);
	
	///Evaluate the interpolation at the given coordinates
	double operator()(double x, double y, double z) const;
private:
	marray<double,1> xcoords;
	marray<double,1> ycoords;
	marray<double,1> zcoords;
	marray<double,3> data;
	marray<double,3> dfdx;
	marray<double,3> dfdy;
	marray<double,3> dfdz;
	marray<double,3> d2fdxdy;
	marray<double,3> d2fdxdz;
	marray<double,3> d2fdydz;
	marray<double,3> d3fdxdydz;
	
	void collect_data(double d[64], std::size_t i, std::size_t j, std::size_t k, double dx, double dy, double dz) const;
protected:
	static void get_coeffs(double c[64], double d[64]);
	static double evaluate_polynomial(double c[64], double t, double u, double v);
};

struct H5File{
    hid_t id;
    H5File(hid_t id):id(id){}
    H5File(const H5File&)=delete;
    H5File(H5File&& h):id(h.id){ h.id=0; }
    H5File& operator=(const H5File&)=delete;
    H5File& operator=(H5File&& h){ std::swap(id,h.id); return *this;}
    ~H5File(){ H5Fclose(id); }
    operator hid_t() const{ return(id); }
};

///\brief A general RAII wrapper around an hid_t
struct H5Handle{
    H5Handle():id(-1),deleter(nullptr){}
    ///\brief Take ownership of an hid_t
    ///\param id The owned HDF5 object reference
    ///\param deleter The HDF5 API function used to dispose of the reference
    ///\param desc A descriptive string used to form an error message if id is not valid
    H5Handle(hid_t id, herr_t(*deleter)(hid_t), const char* desc):
    id(id),deleter(deleter){
        if(id<0){
            if(desc)
                throw std::runtime_error(std::string("Failed to ")+desc);
            else
                throw std::runtime_error("Failed to get HDF5 object");
        }
    }
    H5Handle(const H5Handle&)=delete;
    H5Handle(H5Handle&& other):
    id(other.id),deleter(other.deleter){
        other.id=-1;
    }
    ~H5Handle(){
        if(id>=0)
        (*deleter)(id);
    }
    H5Handle& operator=(H5Handle&)=delete;
    H5Handle& operator=(H5Handle&& other){
        if(&other!=this){
            std::swap(id,other.id);
            std::swap(deleter,other.deleter);
        }
        return *this;
    }
    hid_t get() const{ return id; }
    operator hid_t() const{ return id; }
private:
    hid_t id;
    herr_t(*deleter)(hid_t);
};
	
template<typename T>
struct h5Datatype{};
//not an exhaustive list; other types can be added if necessary
template<> struct h5Datatype<int>{ static hid_t type(){ return H5T_NATIVE_INT; }};
template<> struct h5Datatype<unsigned int>{ static hid_t type(){ return H5T_NATIVE_UINT; }};
template<> struct h5Datatype<float>{ static hid_t type(){ return H5T_NATIVE_FLOAT; }};
template<> struct h5Datatype<double>{ static hid_t type(){ return H5T_NATIVE_DOUBLE; }};
	
bool h5ObjectExists(hid_t loc_id, const char* name);

template<typename DestType, typename SourceType, unsigned int Rank, typename Alloc>
void writeArrayH5(hid_t loc_id, std::string name, const marray<SourceType,Rank,Alloc>& data, unsigned int compressLevel=0){
    hid_t sourceType=h5Datatype<SourceType>::type();
    hid_t destType=h5Datatype<DestType>::type();
    std::array<hsize_t,Rank> extents;
    std::copy_n(data.get_extents().begin(),Rank,extents.begin());
    hid_t dSpace=H5Screate_simple(Rank, &extents[0], nullptr);
    hid_t plist_id=H5P_DEFAULT;
    if(compressLevel){
        assert(compressLevel<=9);
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        std::array<hsize_t,Rank> chunkDims;
        unsigned long size=sizeof(DestType);
        const unsigned long targetChunkSize=1UL<<16; //64KB seems like a good size
        for(unsigned int i=Rank; i>0; i--){
            if(size<targetChunkSize){
                unsigned int dimi=targetChunkSize/size;
                if(dimi>extents[i-1] || dimi==0)
                    dimi=extents[i-1];
                chunkDims[i-1]=dimi;
                size*=chunkDims[i-1];
            }
            else
                chunkDims[i-1]=1;
        }
        H5Pset_chunk(plist_id, Rank, &chunkDims[0]);
        H5Pset_deflate(plist_id, compressLevel);
    }
    
    hid_t dSet=H5Dcreate(loc_id, name.c_str(),destType,dSpace,H5P_DEFAULT,plist_id,H5P_DEFAULT);
    herr_t err=H5Dwrite(dSet,sourceType,H5S_ALL,H5S_ALL,H5P_DEFAULT,&data.front());
    
    //need to close whether or not there's an error
    if(compressLevel)
        H5Pclose(plist_id);
    H5Dclose(dSet);
    H5Sclose(dSpace);
    
    if(err<0)
        throw std::runtime_error("Failed to write marray to HDF5 (name: "
                                 +name+" error: "+std::to_string(err)+")");
}

template<typename T, unsigned int Rank, typename Alloc>
void readArrayH5(hid_t container, const std::string& name, marray<T,Rank,Alloc>& data){
    hid_t array_id=H5Dopen2(container, name.c_str(), H5P_DEFAULT);
    if(array_id<0)
        throw std::runtime_error("Failed to open dataset '"+name+"'");
    hid_t dType=h5Datatype<T>::type();
    hid_t dSpace = H5Dget_space(array_id);
    int raw_dim = H5Sget_simple_extent_ndims(dSpace);
    if(raw_dim!=Rank)
        throw std::runtime_error("Failed to read marray from HDF5: source array has rank "
                                 +std::to_string(raw_dim)+" but target array has rank "
                                 +std::to_string(Rank));
    std::array<hsize_t,Rank> extents;
    H5Sget_simple_extent_dims(dSpace,&extents[0],nullptr);
    H5Sclose(dSpace);
    data.resize(extents);
    herr_t err=H5Dread(array_id,dType,H5S_ALL,H5S_ALL,H5P_DEFAULT,&data.front());
    H5Dclose(array_id);
    if(err<0)
        throw std::runtime_error("Failed to array data from HDF5");
}

template<typename T>
void addH5Attribute(hid_t object, std::string name, const T& contents){
    hid_t dtype=h5Datatype<T>::type();
    hsize_t dim=1;
    hid_t dataspace_id=H5Screate_simple(1, &dim, nullptr);
    hid_t attribute_id=H5Acreate(object,name.c_str(),dtype,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
    H5Awrite(attribute_id,dtype,&contents);
    H5Aclose(attribute_id);
    H5Sclose(dataspace_id);
};

template<>
void addH5Attribute<std::string>(hid_t object, std::string name, const std::string& contents);

template<typename T>
void readH5Attribute(hid_t object, std::string name, T& dest){
    hid_t attribute_id=H5Aopen(object, name.c_str(), H5P_DEFAULT);
    hid_t actualType=H5Aget_type(attribute_id);
    hid_t expectedType=h5Datatype<T>::type();
    
    //TODO: this may be too harsh; it ignores whether conversion is possible
    if(H5Tequal(actualType, expectedType)<=0){
        H5Aclose(attribute_id);
        H5Tclose(actualType);
        throw std::runtime_error("Expected and actual data types for attribute '"+name+"' do not match");
    }
    
    if(H5Aread(attribute_id, expectedType, &dest)<0){
        H5Aclose(attribute_id);
        H5Tclose(actualType);
        throw std::runtime_error("Failed to read attribute '"+name+"'");
    }
    
    H5Aclose(attribute_id);
    H5Tclose(actualType);
}

template<>
void readH5Attribute<std::string>(hid_t object, std::string name, std::string& dest);

} // close namespace

#endif
