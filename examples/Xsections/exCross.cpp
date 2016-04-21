#include "exCross.h"


namespace nusquids{


double NeutrinoDISCrossSectionsFromTablesExtended::LinInter(double x,double xM, double xP,double yM,double yP) const{
  return yM + (yP-yM)*(x-xM)/(xP-xM);
}

double NeutrinoDISCrossSectionsFromTablesExtended::TotalCrossSection(double Enu, NeutrinoFlavor flavor,
                           NeutrinoType neutype, Current current) const{
  // we assume that sterile neutrinos are trully sterile
  if (not (flavor == electron or flavor == muon or flavor == tau))
    return 0.0;

  if (Enu < Emin)
    return 1e-150;

  if (Enu > Emax )
    throw std::runtime_error("NeutrinoCrossSections::Init: Only DIS cross sections are included. Interpolation re\
quested below "+std::to_string(Emin/GeV)+" GeV or above "+std::to_string(Emax/GeV)+" GeV. E_nu = " + std::to_string(Enu/GeV) + " [GeV].");

  // convert to GeV
  Enu /= GeV;

  double logE = log(Enu);
  return gsl_spline_eval(xs_inter[current][neutype][flavor],logE,xs_acc[current][neutype][flavor]);
}

double NeutrinoDISCrossSectionsFromTablesExtended::SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const{
  // we assume that sterile neutrinos are trully sterile
  if (not (flavor == electron or flavor == muon or flavor == tau))
    return 0.0;

  if (E1 < Emin)
    return 1e-150;
  

  if ( E1 < Emin or E1 > Emax )
    throw std::runtime_error("NeutrinoCrossSections::Init: Only DIS cross sections are included. Interpolation re\
quested below "+std::to_string(Emin/GeV)+" GeV or above "+std::to_string(Emax/GeV)+" GeV. E_nu = " + std::to_string(E1/GeV) + " [GeV].");

  // convert to GeV
  E1 /= GeV;
  E2 /= GeV;

  double logE1 = log(E1);
  double logE2 = log(E2);
  double dlogE = logE_data_range[1]-logE_data_range[0];

  size_t loge_M1 = static_cast<size_t>((logE1-logE_data_range[0])/dlogE);
  size_t loge_M2 = static_cast<size_t>((logE2-logE_data_range[0])/dlogE);

  if ( (loge_M2 > div-1) or (loge_M1 > div-1) or (loge_M1 == loge_M2))
    return 1e-150;

  //std::cout << E1 << " " << E2 << " " << loge_M1 << " " << loge_M2 << " " << div << std::endl;
  double phiMM,phiMP,phiPM,phiPP;
  if (current == CC){
    phiMM = dsde_CC_data[loge_M1][loge_M2][neutype][flavor];
    phiMP = dsde_CC_data[loge_M1][loge_M2+1][neutype][flavor];
    if ( loge_M1 == div-1 ){
      // we are at the boundary, cannot bilinearly interpolate
      return LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],
          phiMM,phiMP);
    }
    phiPM = dsde_CC_data[loge_M1+1][loge_M2][neutype][flavor];
    phiPP = dsde_CC_data[loge_M1+1][loge_M2+1][neutype][flavor];
  } else if (current == NC){
    phiMM = dsde_NC_data[loge_M1][loge_M2][neutype][flavor];
    phiMP = dsde_NC_data[loge_M1][loge_M2+1][neutype][flavor];
    if ( loge_M1 == div-1 ){
      // we are at the boundary, cannot bilinearly interpolate
      return LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],
          phiMM,phiMP);
    }
    phiPM = dsde_NC_data[loge_M1+1][loge_M2][neutype][flavor];
    phiPP = dsde_NC_data[loge_M1+1][loge_M2+1][neutype][flavor];
  } else
    throw std::runtime_error("nuSQUIDS::XSECTIONS::ERROR::Current type unkwown.");

  return LinInter(logE1,logE_data_range[loge_M1],logE_data_range[loge_M1+1],
           LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],phiMM,phiMP),
           LinInter(logE2,logE_data_range[loge_M2],logE_data_range[loge_M2+1],phiPM,phiPP));
}

void NeutrinoDISCrossSectionsFromTablesExtended::Init(){
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
          marray<double,2> dsde_CC_raw_data = quickread(filename_dsde_CC);
          marray<double,2> dsde_NC_raw_data = quickread(filename_dsde_NC);
          marray<double,2> sigma_CC_raw_data = quickread(filename_sigma_CC);
          marray<double,2> sigma_NC_raw_data = quickread(filename_sigma_NC);

          // check table shapes and get the number of energy nodes
          unsigned int data_e_size = 0;
          if( sigma_CC_raw_data.extent(0) == sigma_NC_raw_data.extent(0) and
              sigma_NC_raw_data.extent(1) == sigma_NC_raw_data.extent(1) )
            data_e_size = sigma_CC_raw_data.extent(0);
          else
            throw std::runtime_error("nuSQUIDS::xsections::init: Data tables not the same size.");

          // getting the raw data energy node values
          logE_data_range.resize(data_e_size);
          for( int ie = 0; ie < data_e_size; ie ++){
            logE_data_range[ie] = log(sigma_CC_raw_data[ie][0]);
          }

          Emin = sigma_CC_raw_data[0][0]*GeV;
          Emax = sigma_CC_raw_data[data_e_size-1][0]*GeV;
          div = data_e_size;


          // allocate all gsl interpolators
          xs_inter.resize(std::vector<size_t>{2,2,3});
          for ( auto it = xs_inter.begin(); it != xs_inter.end(); it++){
            *it = gsl_spline_alloc(gsl_interp_linear,data_e_size);
          }
          // allocate all gsl interpolators accelerators
          xs_acc.resize(std::vector<size_t>{2,2,3});
          for ( auto it = xs_acc.begin(); it != xs_acc.end(); it++){
            *it = gsl_interp_accel_alloc();
          }

          // initialize gsl interpolators
          // lets loop over active flavors, neutrinos/antineutrinos, and currents.
          for ( Current current : std::vector<Current>{CC,NC}){
            for ( NeutrinoType neutype : std::vector<NeutrinoType>{neutrino,antineutrino}){
              for ( NeutrinoFlavor flavor : std::vector<NeutrinoFlavor>{electron,muon,tau}){
                double sig_data[data_e_size];
                for( unsigned int ie = 0; ie < data_e_size; ie ++){
                  if ( current == CC )
                    sig_data[ie] = sigma_CC_raw_data[ie][1+2*((int)flavor)+(int)neutype];
                  else
                    sig_data[ie] = sigma_NC_raw_data[ie][1+2*((int)flavor)+(int)neutype];
                }
                gsl_spline_init(xs_inter[current][neutype][flavor],logE_data_range.data(),sig_data,data_e_size);
              }
            }
          }

          // convert raw data tables into formatted marrays
          dsde_CC_data.resize(std::vector<size_t>{data_e_size,data_e_size,2,3});
          dsde_NC_data.resize(std::vector<size_t>{data_e_size,data_e_size,2,3});
          for (unsigned int e1 = 0; e1 < data_e_size; e1++){
            for (unsigned int e2 = 0; e2 < data_e_size; e2++){
              for ( NeutrinoType neutype : std::vector<NeutrinoType>{neutrino,antineutrino}){
                for ( NeutrinoFlavor flavor : std::vector<NeutrinoFlavor>{electron,muon,tau}){
                  dsde_CC_data[e1][e2][neutype][flavor] = dsde_CC_raw_data[e1*data_e_size+e2][2+2*(static_cast<int>(flavor))+static_cast<int>(neutype)];
                  dsde_NC_data[e1][e2][neutype][flavor] = dsde_NC_raw_data[e1*data_e_size+e2][2+2*(static_cast<int>(flavor))+static_cast<int>(neutype)];
                }
              }
            }
          }

          /*
          // fill in interpolated cross section tables at the user provided nodes
          marray<double,1> E_range_GeV = logspace(Emin/1.0e9,Emax/1.0e9,div);
          unsigned int e_size = E_range_GeV.size();
          // initialize differential cross section arrays
          dsde_CC_tbl.resize(std::vector<size_t>{e_size,e_size,2,3});
          dsde_NC_tbl.resize(std::vector<size_t>{e_size,e_size,2,3});
          // initialize total cross section arrays
          sigma_CC_tbl.resize(std::vector<size_t>{e_size,2,3});
          sigma_NC_tbl.resize(std::vector<size_t>{e_size,2,3});

          for (unsigned int e1 = 0 ; e1 < e_size ; e1 ++){
              double Enu1 = E_range_GeV[e1];
              // differential cross section
              for (int e2 = 0 ; e2 < e_size ; e2 ++){
                  double Enu2 = E_range_GeV[e2];
                  for ( NeutrinoType neutype : std::vector<NeutrinoType>{neutrino,antineutrino}){
                    for ( NeutrinoFlavor flavor : std::vector<NeutrinoFlavor>{electron,muon,tau}){
                      if ( e2 > e1 ) {
                        dsde_NC_tbl[e1][e2][neutype][flavor] = 0.0;
                        dsde_CC_tbl[e1][e2][neutype][flavor] = 0.0;
                      } else {
                        dsde_NC_tbl[e1][e2][neutype][flavor] = DifferentialCrossSection(Enu1,Enu2,flavor,neutype,NC);
                        dsde_CC_tbl[e1][e2][neutype][flavor] = DifferentialCrossSection(Enu1,Enu2,flavor,neutype,CC);
                      }
                    }
                  }
              }
                  // total cross section

              for ( NeutrinoType neutype : std::vector<NeutrinoType>{neutrino,antineutrino}){
                for ( NeutrinoFlavor flavor : std::vector<NeutrinoFlavor>{electron,muon,tau}){
                  sigma_CC_tbl[e1][neutype][flavor] = TotalCrossSection(Enu1,flavor,neutype,CC);
                  sigma_NC_tbl[e1][neutype][flavor] = TotalCrossSection(Enu1,flavor,neutype,NC);
                }
              }
          }
          */
  } else {
    throw std::runtime_error("nuSQUIDS::XSECTIONS::ERROR::Cross section files not found.");
  }
  // declare the object as initialized;
  is_init = true;
}

NeutrinoDISCrossSectionsFromTablesExtended::~NeutrinoDISCrossSectionsFromTablesExtended(){
  if(is_init){
    // allocate all gsl interpolators
    for ( auto it = xs_inter.begin(); it != xs_inter.end(); it++){
      gsl_spline_free(*it);
    }
    // allocate all gsl interpolators accelerators
    for ( auto it = xs_acc.begin(); it != xs_acc.end(); it++){
      gsl_interp_accel_free(*it);
    }
  }
}


}

