#include "tools.h"

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


bool fexists(std::string filepath)
{
  std::ifstream ifile(filepath.c_str());
  return (bool) ifile;
};

Table quickread(std::string filepath){
    // create and open file stream
    std::ifstream infile(filepath.c_str());

    if(!infile){
        throw std::runtime_error("Error: file could not be opened. Filepath " + filepath);
    }
    #ifdef quickread_DEBUG
    int x,y;
    #endif
    Table table;

    std::string line;
    while(getline(infile,line)){
        Row row;
        std::stringstream linestream(line);

        double data;
        while(linestream >> data){
            row.push_back(data);
        }
        if (!row.empty()){
            #ifdef quickread_DEBUG
            y = row.size();
            #endif
            table.push_back(row);
        }
    }

    #ifdef quickread_DEBUG
    x = table.size();
    cout << "x: " << x << " y: "<< y << endl;
    cout << table[10][0] << endl;
    #endif
    return table;
};

int quickread(std::string filepath,Table tbl){
    // create and open file stream
    std::ifstream infile(filepath.c_str());

    if(!infile){
        std::cerr << "Error: file could not be opened. Filepath " << filepath.c_str()<< std::endl;
        exit(1);
    }
    #ifdef quickread_DEBUG
    int x,y;
    #endif

    std::string line;
    while(getline(infile,line)){
        Row row;
        std::stringstream linestream(line);

        double data;
        while(linestream >> data){
            row.push_back(data);
        }
        if (!row.empty()){
            #ifdef quickread_DEBUG
            y = row.size();
            #endif
            tbl.push_back(row);
        }
    }

    #ifdef quickread_DEBUG
    x = tbl.size();
    cout << "x: " << x << " y: "<< y << endl;
    cout << tbl[10][0] << endl;
    #endif
    return 0;
};

int quickwrite(std::string filepath, Table tbl){
    // create and open file stream
    std::ofstream outfile(filepath.c_str());

    if(!outfile){
        std::cerr << "Error: file could not be created. Filepath " << filepath.c_str()<< std::endl;
        exit(1);
    }

    int tbl_size = tbl.size();
    outfile.precision(15);
    for (int i=0; i < tbl_size;i++){
        std::vector<double> line = tbl[i];
        int line_size = line.size();
        for(int j=0; j < line_size; j++){
           outfile << line[j] << " ";
        }
        outfile << std::endl;
    }

    outfile.close();

    return 0;
};

std::vector<double> linspace(double Emin,double Emax,int div){
    std::vector<double> linpoints;
    double step_lin = (Emax - Emin)/double(div);

    double EE = Emin;
    while (EE <= Emax+0.001){
        linpoints.push_back(EE);
        EE = EE + step_lin;
    }
    return linpoints;
};

std::vector<double> logspace(double Emin,double Emax,int div){
    std::vector<double> logpoints;
    double Emin_log,Emax_log;
    if (Emin < 1.0e-5 ) {
        Emin_log = 0.0;
    } else {
        Emin_log = log(Emin);
    }
    Emax_log = log(Emax);

    double step_log = (Emax_log - Emin_log)/double(div);

    double EE = Emin_log;
    while (EE <= Emax_log+0.001){
        logpoints.push_back(exp(EE));
        EE = EE + step_log;
    }
    return logpoints;
};

Table intertable(Table xy_data, std::vector<double> x_array, int j1 = 0, int j2 = 1){
    Table result;
    int arraysize = xy_data.size();

    double xx[arraysize];
    double yy[arraysize];

    for (int i=0; i < arraysize;i++){
        xx[i] = xy_data[i][j1];
        yy[i] = xy_data[i][j2];
    }

    gsl_spline* inter = gsl_spline_alloc(gsl_interp_cspline,arraysize);
    gsl_interp_accel* inter_accel = gsl_interp_accel_alloc ();
    gsl_spline_init(inter,xx,yy,arraysize);

    for(unsigned int i=0; i < x_array.size();i++){
        Row row;
        row.push_back(x_array[i]);
        if (x_array[i] > xx[arraysize-1] or x_array[i] < xx[0]){
            row.push_back(0.0);
        } else {
            row.push_back(gsl_spline_eval(inter,x_array[i],inter_accel));
        }
        result.push_back(row);
    }

    gsl_spline_free(inter);
    gsl_interp_accel_free(inter_accel);

    return result;
};

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
};

void gsl_matrix_complex_print(gsl_matrix_complex* matrix){
    for(unsigned int i = 0; i < matrix->size1; i++){
        for(unsigned int j = 0; j < matrix->size2; j++){
            std::cout << gsl_matrix_complex_get(matrix,i,j).dat[0] <<
            "+i" << gsl_matrix_complex_get(matrix,i,j).dat[1] << " ";
        }
        std::cout << std::endl;
    }
};

void gsl_matrix_complex_change_basis_UMUC(gsl_matrix_complex* U, gsl_matrix_complex* M){
    int numneu = U->size1;
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
};

void gsl_matrix_complex_change_basis_UCMU(gsl_matrix_complex* U, gsl_matrix_complex* M){
    int numneu = U->size1;
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
};
