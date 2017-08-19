#include <iostream>
#include <iomanip>
#include <vector>
#include <unistd.h>

#include "H5Epublic.h"
#include "H5Tpublic.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Gpublic.h"
#include "H5Fpublic.h"

#include <nuSQuIDS/body.h>

#define H5Gopen_vers 2
#define H5Gcreate_vers 2
#define H5Eset_auto_vers 2

using namespace nusquids;

int main(){
  squids::Const units;

  // open HDF5 file
  hid_t file_id = H5Fcreate("body_test.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t body_group_id,track_group_id;

  // ************************************
  // testing vacuum serialization
  // ************************************
  Vacuum v;
  Vacuum::Track vt(10.*units.km,50*units.km,100.0*units.km);

  // open groups
  hid_t vacuum_group_id = H5Gcreate(file_id, "vacuum", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body_group_id = H5Gcreate(vacuum_group_id, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track_group_id = H5Gcreate(vacuum_group_id, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  v.Serialize(body_group_id);
  vt.Serialize(track_group_id);
  auto vr = Vacuum::Deserialize(body_group_id);
  auto vtr = Vacuum::Track::Deserialize(track_group_id);

  if ( fabs(vr->density(*vtr) - v.density(vt)) >1.0e-5 )
    std::cout << "densities are different after serializing for vacuum" << std::endl;
  if ( fabs(vr->ye(*vtr) - v.ye(vt)) >1.0e-5 )
    std::cout << "ye are different after serializing for vacuum" << std::endl;

  // close hdf5 groups
  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(vacuum_group_id);

  // ************************************
  // testing constant density serialization
  // ************************************
  ConstantDensity c(3.0,0.5);
  ConstantDensity::Track ct(10.*units.km,50*units.km,100.0*units.km);

  // open groups
  hid_t constant_density_group_id = H5Gcreate(file_id, "constant_density", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body_group_id = H5Gcreate(constant_density_group_id, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track_group_id = H5Gcreate(constant_density_group_id, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  c.Serialize(body_group_id);
  ct.Serialize(track_group_id);
  auto cr = ConstantDensity::Deserialize(body_group_id);
  auto ctr = ConstantDensity::Track::Deserialize(track_group_id);

  if ( fabs(cr->density(*ctr) - c.density(ct)) >1.0e-5 )
    std::cout << "densities are different after serializing for constant_density" << std::endl;
  if ( fabs(cr->ye(*ctr) - c.ye(ct)) >1.0e-5 )
    std::cout << "ye are different after serializing for constant_density" << std::endl;

  // close hdf5 groups
  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(constant_density_group_id);

  // ************************************
  // testing variable density serialization
  // ************************************
  std::vector<double> xx {1.,2.,3.,4.,5.};
  std::vector<double> rho {0.,1.,0.,1.,0.};
  std::vector<double> ye {0.5,0.5,0.5,0.5,0.5};
  VariableDensity var(xx,rho,ye);
  VariableDensity::Track vart(10.*units.km,50*units.km,100.0*units.km);

  // open groups
  hid_t variable_density_group_id = H5Gcreate(file_id, "variable_density", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body_group_id = H5Gcreate(variable_density_group_id, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track_group_id = H5Gcreate(variable_density_group_id, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  var.Serialize(body_group_id);
  vart.Serialize(track_group_id);
  auto varr = VariableDensity::Deserialize(body_group_id);
  auto vartr = VariableDensity::Track::Deserialize(track_group_id);

  if ( fabs(varr->density(*ctr) - var.density(vart)) >1.0e-5 )
    std::cout << "densities are different after serializing for constant_density" << std::endl;
  if ( fabs(varr->ye(*vartr) - var.ye(vart)) >1.0e-5 )
    std::cout << "ye are different after serializing for constant_density" << std::endl;

  // close hdf5 groups
  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(variable_density_group_id);

  // ************************************
  // testing earth
  // ************************************
  Earth earth;
  Earth::Track eartht(10.*units.km,50*units.km,100.0*units.km);

  // open groups
  hid_t earth_group_id = H5Gcreate(file_id, "earth", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body_group_id = H5Gcreate(earth_group_id, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track_group_id = H5Gcreate(earth_group_id, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  earth.Serialize(body_group_id);
  eartht.Serialize(track_group_id);
  auto earthr = Earth::Deserialize(body_group_id);
  auto earthtr = Earth::Track::Deserialize(track_group_id);

  if ( fabs(earthr->density(*earthtr) - earth.density(eartht)) >1.0e-5 )
    std::cout << "densities are different after serializing for constant_density" << std::endl;
  if ( fabs(earthr->ye(*earthtr) - earth.ye(eartht)) >1.0e-5 )
    std::cout << "ye are different after serializing for constant_density" << std::endl;

  // close hdf5 groups
  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(earth_group_id);

  // ************************************
  // testing sun
  // ************************************
  Sun sun;
  Sun::Track sunt(10.*units.km,50*units.km,100.0*units.km);

  // open groups
  hid_t sun_group_id = H5Gcreate(file_id, "sun", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body_group_id = H5Gcreate(sun_group_id, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track_group_id = H5Gcreate(sun_group_id, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  sun.Serialize(body_group_id);
  sunt.Serialize(track_group_id);
  auto sunr = Sun::Deserialize(body_group_id);
  auto suntr = Sun::Track::Deserialize(track_group_id);

  if ( fabs(sunr->density(*suntr) - sun.density(sunt)) >1.0e-5 )
    std::cout << "densities are different after serializing for constant_density" << std::endl;
  if ( fabs(sunr->ye(*suntr) - sun.ye(sunt)) >1.0e-5 )
    std::cout << "ye are different after serializing for constant_density" << std::endl;

  // close hdf5 groups
  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(sun_group_id);

  // ************************************
  // testing sunasnu
  // ************************************
  SunASnu sunasnu;
  SunASnu::Track sunasnut(10.*units.km,50*units.km,100.0*units.km);

  // open groups
  hid_t sunasnu_group_id = H5Gcreate(file_id, "sunasnu", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body_group_id = H5Gcreate(sunasnu_group_id, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track_group_id = H5Gcreate(sunasnu_group_id, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  sunasnu.Serialize(body_group_id);
  sunasnut.Serialize(track_group_id);
  auto sunasnur = SunASnu::Deserialize(body_group_id);
  auto sunasnutr = SunASnu::Track::Deserialize(track_group_id);

  if ( fabs(sunasnur->density(*sunasnutr) - sunasnu.density(sunasnut)) >1.0e-5 )
    std::cout << "densities are different after serializing for constant_density" << std::endl;
  if ( fabs(sunasnur->ye(*sunasnutr) - sunasnu.ye(sunasnut)) >1.0e-5 )
    std::cout << "ye are different after serializing for constant_density" << std::endl;

  // close hdf5 groups
  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(sunasnu_group_id);

  // ************************************
  // testing earthatm
  // ************************************
  EarthAtm earthatm;
  EarthAtm::Track earthatmt(10.*units.km,acos(3.1415));

  // open groups
  hid_t earthatm_group_id = H5Gcreate(file_id, "earthatm", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body_group_id = H5Gcreate(earthatm_group_id, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track_group_id = H5Gcreate(earthatm_group_id, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  earthatm.Serialize(body_group_id);
  earthatmt.Serialize(track_group_id);
  auto earthatmr = EarthAtm::Deserialize(body_group_id);
  auto earthatmtr = EarthAtm::Track::Deserialize(track_group_id);

  if ( fabs(earthatmr->density(*earthatmtr) - earthatm.density(earthatmt)) >1.0e-5 )
    std::cout << "densities are different after serializing for constant_density" << std::endl;
  if ( fabs(earthatmr->ye(*earthatmtr) - earthatm.ye(earthatmt)) >1.0e-5 )
    std::cout << "ye are different after serializing for constant_density" << std::endl;

  // close hdf5 groups
  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(earthatm_group_id);

  // closing file
  H5Fclose(file_id);
  H5close();


  return 0;
}
