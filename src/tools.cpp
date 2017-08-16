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

#include "tools.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
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
AkimaSpline::AkimaSpline(const std::vector<double>& x, const std::vector<double>& y){
    assert(x.size()==y.size());
    if(x.size()<3)
        throw std::runtime_error("At least 3 points are required to consturct the Akima spline interpolation");
    const unsigned int n=x.size();
    
    //a lambda, basically, but needs to be able to call itself recursively
    struct mHelper{
        const std::vector<double>& x;
        const std::vector<double>& y;
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
        
        //We mostly resue m values from previous iterations, but the first
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
}

} // close namespace
