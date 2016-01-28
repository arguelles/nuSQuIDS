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

    if(!outfile){
        std::cerr << "Error: file could not be created. Filepath " << filepath.c_str()<< std::endl;
        exit(1);
    }

    outfile.precision(15);
    for (unsigned int i=0; i < tbl.extent(0); i++){
      for ( unsigned int j=0; j < tbl.extent(1); j++){
         outfile << tbl[i][j] << " ";
      }
      outfile << std::endl;
    }

    outfile.close();

    return 0;
}

marray<double,1> linspace(double Emin,double Emax,unsigned int div){
    marray<double,1> linpoints{div+1};
    double step_lin = (Emax - Emin)/double(div);

    double EE = Emin;
    for(unsigned int i=0; i<div; i++, EE+=step_lin)
        linpoints[i] = EE;
    linpoints[div] = Emax;
	
    return linpoints;
}

marray<double,1> logspace(double Emin,double Emax,unsigned int div){
    marray<double,1> logpoints{div+1};
    double Emin_log,Emax_log;
    Emin_log = log(Emin);
    Emax_log = log(Emax);

    double step_log = (Emax_log - Emin_log)/double(div);

	logpoints[0]=Emin;
	double EE = Emin_log+step_log;
	for(unsigned int i=1; i<div; i++, EE+=step_log)
        logpoints[i] = exp(EE);
	logpoints[div]=Emax;
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

} // close namespace
