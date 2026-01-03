#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <unistd.h>

#include "H5Epublic.h"
#include "H5Tpublic.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Gpublic.h"
#include "H5Fpublic.h"

#include <SQuIDS/const.h>
#include <nuSQuIDS/body.h>
#include <nuSQuIDS/resources.h>

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
    std::cout << "densities are different after serializing for constant_density: Before = " << c.density(ct) << " After = " << cr->density(*ctr) << std::endl;
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

  if ( fabs(varr->density(*vartr) - var.density(vart)) >1.0e-5 )
    std::cout << "densities are different after serializing for variable_density: Before = " << var.density(vart) << " After = " << cr->density(*ctr) << std::endl;
  if ( fabs(varr->ye(*vartr) - var.ye(vart)) >1.0e-5 )
    std::cout << "ye are different after serializing for variable_density" << std::endl;

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
    std::cout << "densities are different after serializing for Earth: Before = " << earth.density(eartht) << " After = " << earthr->density(*earthtr) << std::endl;
  if ( fabs(earthr->ye(*earthtr) - earth.ye(eartht)) >1.0e-5 )
    std::cout << "ye are different after serializing for Earth" << std::endl;

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
    std::cout << "densities are different after serializing for Sun: Before = " << sun.density(sunt) << " After = " << sunr->density(*suntr) << std::endl;
  if ( fabs(sunr->ye(*suntr) - sun.ye(sunt)) >1.0e-5 )
    std::cout << "ye are different after serializing for Sun" << std::endl;

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
    std::cout << "densities are different after serializing for SunASnu: Before = " << sunasnu.density(sunasnut) << " After = " << sunasnur->density(*sunasnutr) << std::endl;
  if ( fabs(sunasnur->ye(*sunasnutr) - sunasnu.ye(sunasnut)) >1.0e-5 )
    std::cout << "ye are different after serializing for SunASnu" << std::endl;

  // close hdf5 groups
  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(sunasnu_group_id);

  // ************************************
  // testing earthatm
  // ************************************
  EarthAtm earthatm;
  earthatm.SetAtmosphereHeight(50/*km*/);
  EarthAtm::Track earthatmt(10.*units.km,acos(-1.0),earthatm.GetRadius()*units.km,earthatm.GetAtmosphereHeight()*units.km);

  // open groups
  hid_t earthatm_group_id = H5Gcreate(file_id, "earthatm", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body_group_id = H5Gcreate(earthatm_group_id, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track_group_id = H5Gcreate(earthatm_group_id, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  earthatm.Serialize(body_group_id);
  earthatmt.Serialize(track_group_id);
  auto earthatmr = EarthAtm::Deserialize(body_group_id);
  auto earthatmtr = EarthAtm::Track::Deserialize(track_group_id);

  if (earthatmtr->GetX() != earthatmt.GetX())
    std::cout << "track positions are different after serializing for EarthAtm: Before = " << earthatmt.GetX() << " After = " << earthatmtr->GetX() << std::endl;
  if (earthatmtr->GetInitialX() != earthatmt.GetInitialX())
    std::cout << "track initial positions are different after serializing for EarthAtm: Before = " << earthatmt.GetInitialX() << " After = " << earthatmtr->GetInitialX() << std::endl;
  if (earthatmtr->GetFinalX() != earthatmt.GetFinalX())
    std::cout << "track final positions are different after serializing for EarthAtm: Before = " << earthatmt.GetFinalX() << " After = " << earthatmtr->GetFinalX() << std::endl;
  if ( fabs(earthatmr->density(*earthatmtr) - earthatm.density(earthatmt)) >1.0e-5 )
    std::cout << "densities are different after serializing for EarthAtm: Before = " << earthatm.density(earthatmt) << " After = " << earthatmr->density(*earthatmtr) << std::endl;
  if ( fabs(earthatmr->ye(*earthatmtr) - earthatm.ye(earthatmt)) >1.0e-5 )
    std::cout << "ye are different after serializing for EarthAtm" << std::endl;

  // close hdf5 groups
  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(earthatm_group_id);

  // ************************************
  // testing constant density with composition
  // ************************************
  std::map<PDGCode, double> water_comp;
  water_comp[hydrogen] = 2.0/3.0;
  water_comp[oxygen] = 1.0/3.0;
  ConstantDensity c_comp(1.0, 0.555, water_comp);
  ConstantDensity::Track ct_comp(10.*units.km, 50*units.km, 100.0*units.km);

  hid_t const_comp_group_id = H5Gcreate(file_id, "constant_density_comp", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body_group_id = H5Gcreate(const_comp_group_id, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track_group_id = H5Gcreate(const_comp_group_id, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  c_comp.Serialize(body_group_id);
  ct_comp.Serialize(track_group_id);
  auto cr_comp = ConstantDensity::Deserialize(body_group_id);
  auto ctr_comp = ConstantDensity::Track::Deserialize(track_group_id);

  if (fabs(cr_comp->density(*ctr_comp) - c_comp.density(ct_comp)) > 1.0e-5)
    std::cout << "densities are different after serializing for constant_density with composition" << std::endl;
  if (fabs(cr_comp->ye(*ctr_comp) - c_comp.ye(ct_comp)) > 1.0e-5)
    std::cout << "ye are different after serializing for constant_density with composition" << std::endl;

  // Check composition
  auto orig_comp = c_comp.composition(ct_comp);
  auto read_comp = cr_comp->composition(*ctr_comp);
  if (orig_comp.size() != read_comp.size())
    std::cout << "composition size different after serializing for constant_density: Before = " << orig_comp.size() << " After = " << read_comp.size() << std::endl;
  for (const auto& p : orig_comp) {
    if (read_comp.find(p.first) == read_comp.end())
      std::cout << "composition missing element " << static_cast<int32_t>(p.first) << " after serializing" << std::endl;
    else if (fabs(read_comp.at(p.first) - p.second) > 1.0e-5)
      std::cout << "composition element " << static_cast<int32_t>(p.first) << " different: Before = " << p.second << " After = " << read_comp.at(p.first) << std::endl;
  }

  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(const_comp_group_id);

  // ************************************
  // testing variable density with composition
  // ************************************
  std::vector<double> xx_comp {1., 2., 3., 4., 5.};
  std::vector<double> rho_comp {1., 1., 1., 1., 1.};
  std::vector<double> ye_comp {0.5, 0.5, 0.5, 0.5, 0.5};
  std::map<PDGCode, std::vector<double>> var_composition;
  var_composition[iron] = {0.9, 0.8, 0.7, 0.6, 0.5};
  var_composition[oxygen] = {0.1, 0.2, 0.3, 0.4, 0.5};
  VariableDensity var_comp(xx_comp, rho_comp, ye_comp, var_composition);
  VariableDensity::Track vart_comp(2.0, 2.0, 4.0);  // Position in cm (matching x_arr units)

  hid_t var_comp_group_id = H5Gcreate(file_id, "variable_density_comp", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body_group_id = H5Gcreate(var_comp_group_id, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track_group_id = H5Gcreate(var_comp_group_id, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  var_comp.Serialize(body_group_id);
  vart_comp.Serialize(track_group_id);
  auto varr_comp = VariableDensity::Deserialize(body_group_id);
  auto vartr_comp = VariableDensity::Track::Deserialize(track_group_id);

  // Check composition at a point
  auto var_orig_comp = var_comp.composition(vart_comp);
  auto var_read_comp = varr_comp->composition(*vartr_comp);
  if (var_orig_comp.size() != var_read_comp.size())
    std::cout << "composition size different after serializing for variable_density: Before = " << var_orig_comp.size() << " After = " << var_read_comp.size() << std::endl;
  for (const auto& p : var_orig_comp) {
    if (var_read_comp.find(p.first) == var_read_comp.end())
      std::cout << "variable_density composition missing element " << static_cast<int32_t>(p.first) << " after serializing" << std::endl;
    else if (fabs(var_read_comp.at(p.first) - p.second) > 1.0e-5)
      std::cout << "variable_density composition element " << static_cast<int32_t>(p.first) << " different: Before = " << p.second << " After = " << var_read_comp.at(p.first) << std::endl;
  }

  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(var_comp_group_id);

  // ************************************
  // testing earthatm with composition
  // ************************************
  std::string prem_path = getResourcePath() + "/astro/EARTH_MODEL_PREM_wIso.dat";
  EarthAtm earthatm_comp(prem_path);
  earthatm_comp.SetAtmosphereHeight(50/*km*/);
  EarthAtm::Track earthatmt_comp = earthatm_comp.MakeTrackWithCosine(-1.0);  // Vertical through Earth
  earthatmt_comp.SetX(0.5 * earthatmt_comp.GetFinalX());  // Midpoint (core)

  hid_t earthatm_comp_group_id = H5Gcreate(file_id, "earthatm_comp", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  body_group_id = H5Gcreate(earthatm_comp_group_id, "body", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  track_group_id = H5Gcreate(earthatm_comp_group_id, "track", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  earthatm_comp.Serialize(body_group_id);
  earthatmt_comp.Serialize(track_group_id);
  auto earthatmr_comp = EarthAtm::Deserialize(body_group_id);
  auto earthatmtr_comp = EarthAtm::Track::Deserialize(track_group_id);

  if (fabs(earthatmr_comp->density(*earthatmtr_comp) - earthatm_comp.density(earthatmt_comp)) > 1.0e-5)
    std::cout << "densities are different after serializing for EarthAtm with composition: Before = " << earthatm_comp.density(earthatmt_comp) << " After = " << earthatmr_comp->density(*earthatmtr_comp) << std::endl;
  if (fabs(earthatmr_comp->ye(*earthatmtr_comp) - earthatm_comp.ye(earthatmt_comp)) > 1.0e-5)
    std::cout << "ye are different after serializing for EarthAtm with composition" << std::endl;

  // Check composition at Earth's core (should have significant iron)
  auto earth_orig_comp = earthatm_comp.composition(earthatmt_comp);
  auto earth_read_comp = earthatmr_comp->composition(*earthatmtr_comp);
  if (earth_orig_comp.size() != earth_read_comp.size())
    std::cout << "composition size different after serializing for EarthAtm: Before = " << earth_orig_comp.size() << " After = " << earth_read_comp.size() << std::endl;
  for (const auto& p : earth_orig_comp) {
    if (earth_read_comp.find(p.first) == earth_read_comp.end())
      std::cout << "EarthAtm composition missing element " << static_cast<int32_t>(p.first) << " after serializing" << std::endl;
    else if (fabs(earth_read_comp.at(p.first) - p.second) > 1.0e-5)
      std::cout << "EarthAtm composition element " << static_cast<int32_t>(p.first) << " different: Before = " << p.second << " After = " << earth_read_comp.at(p.first) << std::endl;
  }

  // Verify iron is present in core
  if (earth_read_comp.find(iron) == earth_read_comp.end() || earth_read_comp.at(iron) < 0.5)
    std::cout << "EarthAtm composition should have significant iron at core after serialization" << std::endl;

  H5Gclose(body_group_id);
  H5Gclose(track_group_id);
  H5Gclose(earthatm_comp_group_id);

  // closing file
  H5Fclose(file_id);
  H5close();


  return 0;
}
