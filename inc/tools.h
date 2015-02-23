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

// file array
typedef std::vector<double> Row;
typedef std::vector<Row> Table;

// string convertion
//std::string toString(double);
//std::string toString(int);

// file read and write
bool fexists(std::string);
Table quickread(std::string);
int quickwrite(std::string,Table&);
/// \brief Construct a linear space
/// @param min Minimum value in the linear span.
/// @param max Maximum value in the linear span.
/// @param div Number of divisions in the span.
marray<double,1> linspace(double min,double max,unsigned int div);
/// \brief 
marray<double,1> logspace(double,double,unsigned int);
Table intertable(Table&,std::vector<double>,unsigned int,unsigned int);
// additional GSL-like tools
void gsl_matrix_complex_conjugate(gsl_matrix_complex*);
void gsl_matrix_complex_print(gsl_matrix_complex*);
void gsl_matrix_complex_change_basis_UMUC(gsl_matrix_complex*, gsl_matrix_complex*);
void gsl_matrix_complex_change_basis_UCMU(gsl_matrix_complex*, gsl_matrix_complex*);

} // close namespace

#endif
