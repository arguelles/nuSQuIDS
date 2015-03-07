#include "tools.h"

namespace nusquids{

/*
std::string toString(double value) {
   std::stringstream ss;
   ss << value;
   return ss.str();
}

std::string toString(int value) {
   std::stringstream ss;
   ss << value;
   return ss.str();
}
*/

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
    unsigned int i = 0;
    while (EE <= Emax+0.001){
        linpoints[i] = EE;
        EE = EE + step_lin;
        i++;
    }
    return linpoints;
}

marray<double,1> logspace(double Emin,double Emax,unsigned int div){
    marray<double,1> logpoints{div+1};
    double Emin_log,Emax_log;
    if (Emin < 1.0e-5 ) {
        Emin_log = 0.0;
    } else {
        Emin_log = log(Emin);
    }
    Emax_log = log(Emax);

    double step_log = (Emax_log - Emin_log)/double(div);

    double EE = Emin_log;
    unsigned int i = 0;
    while (EE <= Emax_log+0.001){
        logpoints[i] = exp(EE);
        EE = EE + step_log;
        i++;
    }
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
