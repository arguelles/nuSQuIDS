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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <stdexcept>
#include "marray.h"

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
/// @param div Number of divisions in the span.
marray<double,1> linspace(double min,double max,unsigned int div);
/// \brief Construct a logarithmic space.
/// @param min Minimum value in the logarithmic span.
/// @param max Maximum value in the logarithmic span.
/// @param div Number of divisions in the span.
marray<double,1> logspace(double min,double max,unsigned int div);
// additional GSL-like tools
void gsl_matrix_complex_conjugate(gsl_matrix_complex*);
void gsl_matrix_complex_print(gsl_matrix_complex*);
void gsl_matrix_complex_change_basis_UMUC(gsl_matrix_complex*, gsl_matrix_complex*);
void gsl_matrix_complex_change_basis_UCMU(gsl_matrix_complex*, gsl_matrix_complex*);

} // close namespace

#endif
