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


#ifndef __TOOLS_H
#define __TOOLS_H

#include "version.h"
#include <algorithm>
#include <string>
#include <vector>
#include "marray.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>

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
// additional GSL-like tools
void gsl_matrix_complex_conjugate(gsl_matrix_complex*);
void gsl_matrix_complex_print(gsl_matrix_complex*);
void gsl_matrix_complex_change_basis_UMUC(gsl_matrix_complex*, gsl_matrix_complex*);
void gsl_matrix_complex_change_basis_UCMU(gsl_matrix_complex*, gsl_matrix_complex*);
	
class AkimaSpline{
private:
    struct segment{
        double x;
        double a0, a1, a2, a3;
        
        double operator()(double x) const{
            x-=this->x;
            return((((a3*x)+a2)*x+a1)*x+a0);
        }
        operator double() const{ return(x); }
    };
    std::vector<segment> segments;
public:
    AkimaSpline(){}
    AkimaSpline(const std::vector<double>& x, const std::vector<double>& y);
    
    double operator()(double x) const{
        if(x<segments.front().x)
            return(segments.front()(x));
        auto seg_it=std::upper_bound(segments.begin(),segments.end(),x);
        seg_it--;
        return((*seg_it)(x));
    }
};

} // close namespace

#endif
