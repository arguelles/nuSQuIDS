#include "xsections.h"

namespace nusquids{

NeutrinoCrossSections::NeutrinoCrossSections(){};

double NeutrinoCrossSections::LinInter(double x,double xM, double xP,double yM,double yP){
  return yM + (yP-yM)*(x-xM)/(xP-xM);
}

double NeutrinoCrossSections::SigmaInter(double logE, gsl_spline * inter, gsl_interp_accel * acc){
  return gsl_spline_eval(inter,logE,acc);
}

double NeutrinoCrossSections::DSDEInter(double Enu1,double Enu2,int flavor,std::vector<double> loge_range, std::string current){
  double logE1 = log(Enu1);
  double logE2 = log(Enu2);
  double dlogE = loge_range[1]-loge_range[0];

  int loge_M1 = (int)((logE1-loge_range[0])/dlogE);
  int loge_M2 = (int)((logE2-loge_range[0])/dlogE);

  if (loge_M2 < 0 or loge_M1 < 0)
    return 0.0;

  double phiMM,phiMP,phiPM,phiPP;
  if (current == "CC"){
    phiMM = dsde_CC_data[loge_M1*loge_range.size()+loge_M2][2+flavor];
    phiMP = dsde_CC_data[loge_M1*loge_range.size()+(loge_M2+1)][2+flavor];
    phiPM = dsde_CC_data[(loge_M1+1)*loge_range.size()+(loge_M2)][2+flavor];
    phiPP = dsde_CC_data[(loge_M1+1)*loge_range.size()+(loge_M2+1)][2+flavor];
  } else if (current == "NC"){
    phiMM = dsde_NC_data[loge_M1*loge_range.size()+loge_M2][2+flavor];
    phiMP = dsde_NC_data[loge_M1*loge_range.size()+(loge_M2+1)][2+flavor];
    phiPM = dsde_NC_data[(loge_M1+1)*loge_range.size()+(loge_M2)][2+flavor];
    phiPP = dsde_NC_data[(loge_M1+1)*loge_range.size()+(loge_M2+1)][2+flavor];
  } else
    throw std::runtime_error("nuSQUIDS::XSECTIONS::ERROR::Current type unkwown.");

  return LinInter(logE1,loge_range[loge_M1],loge_range[loge_M1+1],
           LinInter(logE2,loge_range[loge_M2],loge_range[loge_M2+1],phiMM,phiMP),
           LinInter(logE2,loge_range[loge_M2],loge_range[loge_M2+1],phiPM,phiPP));
}

void NeutrinoCrossSections::Init(double Emin_in, double Emax_in, int div_in){

       Emin = Emin_in;
       Emax = Emax_in;
       div = div_in;

       std::string root = XSECTION_LOCATION ;
       std::string filename_format = "_1e+11_1e+18_500.dat";

       std::string filename_dsde_CC = root+"dsde_CC"+filename_format;
       std::string filename_dsde_NC = root+"dsde_NC"+filename_format;
       std::string filename_sigma_CC = root+"sigma_CC"+filename_format;
       std::string filename_sigma_NC = root+"sigma_NC"+filename_format;

       // check if files exist for this energies and divisions
       if(
          fexists(filename_dsde_CC) and
          fexists(filename_dsde_NC) and
          fexists(filename_sigma_CC) and
          fexists(filename_sigma_NC)
          )
       {
          // read data tables
          dsde_CC_data = quickread(filename_dsde_CC);
          dsde_NC_data = quickread(filename_dsde_NC);
          sigma_CC_data = quickread(filename_sigma_CC);
          sigma_NC_data = quickread(filename_sigma_NC);

          int data_e_size = sigma_CC_data.size();
          std::vector<gsl_spline *> xs_cc_inter(6);
          std::vector<gsl_interp_accel *> xs_cc_inter_accel(6);
          std::vector<gsl_spline *> xs_nc_inter(6);
          std::vector<gsl_interp_accel *> xs_nc_inter_accel(6);

          std::vector<double> logE_data_range(data_e_size);
          for( int ie = 0; ie < data_e_size; ie ++){
            logE_data_range[ie] = log(sigma_CC_data[ie][0]);
          }

          for ( int flavor = 0; flavor < 6; flavor ++ ){
            xs_cc_inter[flavor] = gsl_spline_alloc(gsl_interp_linear,data_e_size);
            xs_cc_inter_accel[flavor] = gsl_interp_accel_alloc();

            xs_nc_inter[flavor] = gsl_spline_alloc(gsl_interp_linear,data_e_size);
            xs_nc_inter_accel[flavor] = gsl_interp_accel_alloc();

            double sig_cc_data[data_e_size], sig_nc_data[data_e_size];
            for( int ie = 0; ie < data_e_size; ie ++){
              sig_cc_data[ie] = sigma_CC_data[ie][1+flavor];
              sig_nc_data[ie] = sigma_NC_data[ie][1+flavor];
            }
            gsl_spline_init(xs_cc_inter[flavor],logE_data_range.data(),sig_cc_data,data_e_size);
            gsl_spline_init(xs_nc_inter[flavor],logE_data_range.data(),sig_nc_data,data_e_size);
          }

          // fill in interpolated cross section tables
          marray<double,1> E_range_GeV = logspace(Emin/1.0e9,Emax/1.0e9,div);
          int e_size = E_range_GeV.size();

          for (int e1 = 0 ; e1 < e_size ; e1 ++){
              double Enu1 = E_range_GeV[e1];
              for (int e2 = 0 ; e2 < e_size ; e2 ++){
                  double Enu2 = E_range_GeV[e2];
                  Row dsde_cc;
                  Row dsde_nc;

                  dsde_cc.push_back(Enu1);
                  dsde_nc.push_back(Enu1);
                  dsde_cc.push_back(Enu2);
                  dsde_nc.push_back(Enu2);

                  for (int flavor = 0 ; flavor < 6; flavor ++){
                      if ( e2 > e1 ) {
                        dsde_cc.push_back(0.0);
                        dsde_nc.push_back(0.0);
                      } else {
                        dsde_cc.push_back(DSDEInter(Enu1,Enu2,flavor,logE_data_range,"CC"));
                        dsde_nc.push_back(DSDEInter(Enu1,Enu2,flavor,logE_data_range,"NC"));
                      }
                  }
                  dsde_CC_tbl.push_back(dsde_cc);
                  dsde_NC_tbl.push_back(dsde_nc);
              }
              Row row_cc;
              Row row_nc;

              row_cc.push_back(Enu1);
              row_nc.push_back(Enu1);

              for (int flavor = 0 ; flavor < 6; flavor ++){
                  row_cc.push_back(SigmaInter(log(Enu1),xs_cc_inter[flavor],xs_cc_inter_accel[flavor]));
                  row_nc.push_back(SigmaInter(log(Enu1),xs_nc_inter[flavor],xs_nc_inter_accel[flavor]));
              }

              sigma_CC_tbl.push_back(row_cc);
              sigma_NC_tbl.push_back(row_nc);
          }

         for(int i = 0; i < 6; i++){
          free(xs_cc_inter[i]);
          free(xs_nc_inter[i]);
          free(xs_cc_inter_accel[i]);
          free(xs_nc_inter_accel[i]);
         }
       } else {
          throw std::runtime_error("nuSQUIDS::XSECTIONS::ERROR::Cross section files not found.");
       }

};

NeutrinoCrossSections::NeutrinoCrossSections(double Emin_in,double Emax_in,int div_in){
  Init(Emin_in,Emax_in,div_in);
};

// differential cross sections

double NeutrinoCrossSections::dsde_CC(int i_enu, int i_ele, int flv, int neutype){
    if (i_ele > i_enu) {
        return 0.0;
    } else {
        int ii = i_enu*(div+1) + (i_ele);
        return dsde_CC_tbl[ii][2+2*flv+neutype];
    }
};

double NeutrinoCrossSections::dsde_NC(int i_enu, int i_ele, int flv, int neutype){
    if (i_ele > i_enu) {
        return 0.0;
    } else {
        int ii = i_enu*(div+1) + i_ele;
        return (double) dsde_NC_tbl[ii][2+neutype];
    }
};

//total cross sections

double NeutrinoCrossSections::sigma_CC(int i_enu, int flv, int neutype){
    return sigma_CC_tbl[i_enu][1+2*flv+neutype];
};

double NeutrinoCrossSections::sigma_NC(int i_enu, int flv, int neutype){
    return sigma_NC_tbl[i_enu][1+2*flv+neutype];
};

} // close namespace
