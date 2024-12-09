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


#include <nuSQuIDS/body.h>

#include <cmath>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <H5Apublic.h>
#include <H5LTpublic.h>
#include <H5Ppublic.h>
#include <H5Spublic.h>
#include <H5public.h>
#include <H5version.h>

#include <SQuIDS/const.h>

#include <nuSQuIDS/resources.h>

// Macros
#define SQR(x)      ((x)*(x))                        // x^2

namespace nusquids{

namespace {
  
squids::Const param;

void addStringAttribute(hid_t object, std::string name, std::string contents){
  hid_t strtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(strtype, contents.size());
  hsize_t dim=1;
  hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
  hid_t attribute_id = H5Acreate(object,name.c_str(),strtype,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
  H5Awrite(attribute_id, strtype, &contents[0]);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);
}

//not currently used
/*
std::string readStringAttribute(hid_t object, std::string name){
  hid_t strtype = H5Tcopy(H5T_C_S1);
  hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
  hsize_t storage = H5Aget_storage_size(attribute_id);
  if(storage==0)
    throw std::runtime_error("Not finite space");
  std::unique_ptr<char[]> char_out(new char[storage]);
  herr_t status = H5Aread(attribute_id,strtype,char_out.get());
  if(status<0)
    throw std::runtime_error("Failed to read attribute '"+name+"'");
  H5Aclose(attribute_id);
  return std::string(char_out.get());
}
*/

void addDoubleAttribute(hid_t object, std::string name, double value){
  hsize_t dim=1;
  hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
  hid_t attribute_id = H5Acreate(object,name.c_str(),H5T_IEEE_F64LE,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_IEEE_F64LE, &value);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);
}

double readDoubleAttribute(hid_t object, std::string name){
  double target;
  hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
  herr_t status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &target);
  if(status<0)
    throw std::runtime_error("Failed to read attribute '"+name+"'");
  H5Aclose(attribute_id);
  return target;
}

void addUIntAttribute(hid_t object, std::string name, unsigned int value){
  hsize_t dim=1;
  hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
  hid_t attribute_id = H5Acreate(object,name.c_str(),H5T_NATIVE_UINT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_UINT, &value);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);
}

unsigned int readUIntAttribute(hid_t object, std::string name){
  unsigned int target;
  hid_t attribute_id = H5Aopen(object,name.c_str(),H5P_DEFAULT);
  herr_t status = H5Aread(attribute_id, H5T_NATIVE_UINT, &target);
  if(status<0)
    throw std::runtime_error("Failed to read attribute '"+name+"'");
  H5Aclose(attribute_id);
  return target;
}

} // close unnamed namespace


/*
----------------------------------------------------------------------
         VACUUM CLASS DEFINITIONS
----------------------------------------------------------------------
*/

double Vacuum::density(const GenericTrack& track_input) const{
  return 0.0;
}

double Vacuum::ye(const GenericTrack& track_input) const{
  return 1.0;
}

bool Vacuum::IsConstantDensity() const { return true;}

void Vacuum::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
}

std::shared_ptr<Vacuum> Vacuum::Deserialize(hid_t group){
  return std::make_shared<Vacuum>();
}

void Vacuum::Track::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addDoubleAttribute(group,"x",x);
  addDoubleAttribute(group,"xini",xini);
  addDoubleAttribute(group,"xend",xend);
}

std::shared_ptr<Vacuum::Track> Vacuum::Track::Deserialize(hid_t group){
  double x_   =readDoubleAttribute(group,"x");
  double xini_=readDoubleAttribute(group,"xini");
  double xend_=readDoubleAttribute(group,"xend");
  return std::make_shared<Vacuum::Track>(x_,xini_,xend_);
}

/*
----------------------------------------------------------------------
         ConstantDensity CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
ConstantDensity::ConstantDensity(double constant_density,double constant_ye):
Body(),
constant_density(constant_density),
constant_ye(constant_ye)
{
  BodyParams = {constant_density, constant_ye};
}

void ConstantDensity::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addDoubleAttribute(group,"constant_density",constant_density);
  addDoubleAttribute(group,"constant_ye",constant_ye);
}

std::shared_ptr<ConstantDensity> ConstantDensity::Deserialize(hid_t group){
  double const_dens=readDoubleAttribute(group,"constant_density");
  double const_ye=readDoubleAttribute(group,"constant_ye");
  return std::make_shared<ConstantDensity>(const_dens,const_ye);
}

// track constructor

void ConstantDensity::Track::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addDoubleAttribute(group,"x",x);
  addDoubleAttribute(group,"xini",xini);
  addDoubleAttribute(group,"xend",xend);
}

std::shared_ptr<ConstantDensity::Track> ConstantDensity::Track::Deserialize(hid_t group){
  double x_   =readDoubleAttribute(group,"x");
  double xini_=readDoubleAttribute(group,"xini");
  double xend_=readDoubleAttribute(group,"xend");
  return std::make_shared<ConstantDensity::Track>(x_,xini_,xend_);
}

double ConstantDensity::density(const GenericTrack& track_input) const
{
  return constant_density;
}

double ConstantDensity::ye(const GenericTrack& track_input) const
{
  return constant_ye;
}

bool ConstantDensity::IsConstantDensity() const { return true;}

/*
----------------------------------------------------------------------
         VariableDensity CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
VariableDensity::VariableDensity(std::vector<double> x_input,std::vector<double> density_input,std::vector<double> ye_input):
Body(),x_arr(std::move(x_input)),density_arr(std::move(density_input)),ye_arr(std::move(ye_input)),
inter_density(x_arr,density_arr),inter_ye(x_arr,ye_arr)
{
  assert("nuSQUIDS::Error::VariableDensityConstructor: Invalid array sizes." && x_arr.size() == density_arr.size() && x_arr.size() == ye_arr.size());
  arraysize = x_arr.size();

  x_min = x_arr.front();
  x_max = x_arr.back();

  for(double xx : x_arr)
    BodyParams.push_back(xx);
  for(double rho : density_arr)
    BodyParams.push_back(rho);
  for(double ye : ye_arr)
    BodyParams.push_back(ye);
}

void VariableDensity::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addUIntAttribute(group,"arraysize",arraysize);
  std::vector<hsize_t> dims {arraysize};
  H5LTmake_dataset_double(group, "x_arr", 1, dims.data(), x_arr.data());
  H5LTmake_dataset_double(group, "density_arr", 1, dims.data(), density_arr.data());
  H5LTmake_dataset_double(group, "ye_arr", 1, dims.data(), ye_arr.data());
}

std::shared_ptr<VariableDensity> VariableDensity::Deserialize(hid_t group){
  unsigned int asize=readUIntAttribute(group,"arraysize");
  std::vector<double> x_vec(asize),rho_vec(asize),ye_vec(asize);
  H5LTread_dataset_double(group,"x_arr",x_vec.data());
  H5LTread_dataset_double(group,"density_arr",rho_vec.data());
  H5LTread_dataset_double(group,"ye_arr",ye_vec.data());
  return std::make_shared<VariableDensity>(x_vec,rho_vec,ye_vec);
}

// track constructor

double VariableDensity::density(const GenericTrack& track_input) const
{
  double x = track_input.GetX()/param.cm;
  if (x < x_min or x > x_max ){
    return 0;
  } else {
    return inter_density(x);
  }
}
double VariableDensity::ye(const GenericTrack& track_input) const
{
  double x = track_input.GetX()/param.cm;
  if (x < x_min or x > x_max ){
    return 0;
  } else {
    return inter_ye(x);
  }
}

void VariableDensity::Track::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addDoubleAttribute(group,"x",x);
  addDoubleAttribute(group,"xini",xini);
  addDoubleAttribute(group,"xend",xend);
}

std::shared_ptr<VariableDensity::Track> VariableDensity::Track::Deserialize(hid_t group){
  double x_   =readDoubleAttribute(group,"x");
  double xini_=readDoubleAttribute(group,"xini");
  double xend_=readDoubleAttribute(group,"xend");
  return std::make_shared<VariableDensity::Track>(x_,xini_,xend_);
}

VariableDensity::~VariableDensity(){}

/*
----------------------------------------------------------------------
         Earth CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
Earth::Earth():Earth(getResourcePath()+"/astro/EARTH_MODEL_PREM.dat"){}

Earth::Earth(std::string filepath):Body()
{
  // The Input file should have the radius specified from 0 to 1.
  // where 0 is the center of the Earth and 1 is the surface.
  radius = 6371.0; // [km]

  marray<double,2> earth_model = quickread(filepath);
  arraysize = earth_model.extent(0);

  earth_radius.resize(arraysize);
  earth_density.resize(arraysize);
  earth_ye.resize(arraysize);

  for (size_t i = 0; i < arraysize; i++){
    earth_radius[i] = earth_model[i][0];
    earth_density[i] = earth_model[i][1];
    earth_ye[i] = earth_model[i][2];
  }

  x_radius_min = earth_radius[0];
  x_radius_max = earth_radius[arraysize-1];
  x_rho_min = earth_density[0];
  x_rho_max = earth_density[arraysize-1];
  x_ye_min = earth_ye[0];
  x_ye_max = earth_ye[arraysize-1];

  inter_density=AkimaSpline(earth_radius,earth_density);
  inter_ye=AkimaSpline(earth_radius,earth_ye);
}

Earth::Earth(std::vector<double> x,std::vector<double> rho,std::vector<double> ye):
Body(),earth_radius(std::move(x)),earth_density(std::move(rho)),earth_ye(std::move(ye)),
inter_density(earth_radius,earth_density),inter_ye(earth_radius,earth_ye)
{
  assert("nuSQUIDS::Error::EarthConstructor: Invalid array sizes." && earth_radius.size() == earth_density.size() && earth_radius.size() == earth_ye.size());
  // The Input file should have the radius specified from 0 to 1.
  // where 0 is the center of the Earth and 1 is the surface.
  radius = 6371.0; // [km]
  arraysize = earth_radius.size();

  x_radius_min = earth_radius[0];
  x_radius_max = earth_radius[arraysize-1];
  x_rho_min = earth_density[0];
  x_rho_max = earth_density[arraysize-1];
  x_ye_min = earth_ye[0];
  x_ye_max = earth_ye[arraysize-1];
}

void Earth::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addUIntAttribute(group,"arraysize",arraysize);
  std::vector<hsize_t> dims {arraysize};
  H5LTmake_dataset_double(group, "earth_radius", 1, dims.data(), earth_radius.data());
  H5LTmake_dataset_double(group, "earth_density", 1, dims.data(), earth_density.data());
  H5LTmake_dataset_double(group, "earth_ye", 1, dims.data(), earth_ye.data());
}

std::shared_ptr<Earth> Earth::Deserialize(hid_t group){
  unsigned int asize=readUIntAttribute(group,"arraysize");
  std::vector<double> x_vec(asize),rho_vec(asize),ye_vec(asize);
  H5LTread_dataset_double(group,"earth_radius",x_vec.data());
  H5LTread_dataset_double(group,"earth_density",rho_vec.data());
  H5LTread_dataset_double(group,"earth_ye",ye_vec.data());
  return std::make_shared<Earth>(x_vec,rho_vec,ye_vec);
}

double Earth::density(const GenericTrack& track_input) const
{
  const Earth::Track& track_earth = static_cast<const Earth::Track&>(track_input);
  double xkm = track_earth.GetX()/param.km;
  double r2 = SQR(radius)+SQR(xkm)-(track_earth.GetBaseline()/param.km)*xkm;
  double r;
  if (r2 > 0.0)
    r = sqrt(r2);
  else if(fabs(r2) < 1.e-6)
    r = 0;
  else
    throw std::runtime_error("nuSQUIDS::Earth::ye got impossible geometry.");

  if ( r/radius < x_radius_min ){
    return x_rho_min;
  }
  else if ( r/radius > x_radius_max ) {
    return x_rho_max;
  }
  else {
    return inter_density(r/radius);
  }
}

double Earth::ye(const GenericTrack& track_input) const
{
  const Earth::Track& track_earth = static_cast<const Earth::Track&>(track_input);
  double xkm = track_earth.GetX()/param.km;
  double r2 = SQR(radius)+SQR(xkm)-(track_earth.GetBaseline()/param.km)*xkm;
  double r;
  if (r2 > 0.0)
    r = sqrt(r2);
  else if(fabs(r2) < 1.e-6)
    r = 0;
  else
    throw std::runtime_error("nuSQUIDS::Earth::ye got impossible geometry.");

  if ( r/radius < x_radius_min ){
    return x_ye_min;
  }
  else if ( r/radius > x_radius_max ) {
    return x_ye_max;
  }
  else {
    return inter_ye(r/radius);
  }
}

Earth::~Earth(){}

// track constructor
void Earth::Track::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addDoubleAttribute(group,"x",x);
  addDoubleAttribute(group,"xini",xini);
  addDoubleAttribute(group,"xend",xend);
  addDoubleAttribute(group,"baseline",baseline);
}

std::shared_ptr<Earth::Track> Earth::Track::Deserialize(hid_t group){
  double x_   =readDoubleAttribute(group,"x");
  double xini_=readDoubleAttribute(group,"xini");
  double xend_=readDoubleAttribute(group,"xend");
  double baseline_=readDoubleAttribute(group,"baseline");
  return std::make_shared<Earth::Track>(x_,xini_,xend_,baseline_);
}

void Earth::Track::FillDerivedParams(std::vector<double>& TrackParams) const{
  TrackParams.push_back(baseline);
}

/*
----------------------------------------------------------------------
         SUN CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
Sun::Sun():Sun(getResourcePath()+"/astro/bs05_agsop.dat")
{}

Sun::Sun(std::string sunlocation):Body()
{
  radius = 695980.0*param.km;

  // import sun model
  sun_model = quickread(sunlocation);
  arraysize = sun_model.extent(0);

  sun_radius.resize(arraysize);
  sun_density.resize(arraysize);
  sun_xh.resize(arraysize);

  for (unsigned int i=0; i < arraysize;i++){
    sun_radius[i] = sun_model[i][1];
    sun_density[i] = sun_model[i][3];
    sun_xh[i] = sun_model[i][6];
  }

  inter_density=AkimaSpline(sun_radius,sun_density);
  inter_xh=AkimaSpline(sun_radius,sun_xh);
}

Sun::Sun(std::vector<double> x,std::vector<double> rho,std::vector<double> xh):
Body(),sun_radius(std::move(x)),sun_density(std::move(rho)),sun_xh(std::move(xh)),
inter_density(sun_radius,sun_density),inter_xh(sun_radius,sun_xh)
{
  radius = 695980.0*param.km;
  arraysize = sun_radius.size();
}

void Sun::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addUIntAttribute(group,"arraysize",arraysize);
  std::vector<hsize_t> dims {arraysize};
  H5LTmake_dataset_double(group, "sun_radius", 1, dims.data(), sun_radius.data());
  H5LTmake_dataset_double(group, "sun_density", 1, dims.data(), sun_density.data());
  H5LTmake_dataset_double(group, "sun_xh", 1, dims.data(), sun_xh.data());
}
std::shared_ptr<Sun> Sun::Deserialize(hid_t group){
  unsigned int asize=readUIntAttribute(group,"arraysize");
  std::vector<double> x_vec(asize),rho_vec(asize),xh_vec(asize);
  H5LTread_dataset_double(group,"sun_radius",x_vec.data());
  H5LTread_dataset_double(group,"sun_density",rho_vec.data());
  H5LTread_dataset_double(group,"sun_xh",xh_vec.data());
  return std::make_shared<Sun>(x_vec,rho_vec,xh_vec);
}

double Sun::rdensity(double x) const{
// x is adimentional radius : x = 0 : center, x = 1 : radius
  if (x < sun_radius[0]){
    return sun_density[0];
  } else if ( x > sun_radius[arraysize-1]){
    return 0;
  } else {
    return inter_density(x);
  }
}

double Sun::rxh(double x) const{
// x is adimentional radius : x = 0 : center, x = 1 : radius
  if (x < sun_radius[0]){
    return sun_xh[0];
  } else if ( x > sun_radius[arraysize-1]){
    return 0;
  } else {
    return inter_xh(x);
  }
}

double Sun::density(const GenericTrack& track_input) const
{
  double r = track_input.GetX()/(radius);
  return rdensity(r);
}
double Sun::ye(const GenericTrack& track_input) const
{
  double r = track_input.GetX()/(radius);
  return 0.5*(1.0+rxh(r));
}

Sun::~Sun(){}

void Sun::Track::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addDoubleAttribute(group,"x",x);
  addDoubleAttribute(group,"xini",xini);
  addDoubleAttribute(group,"xend",xend);
}

std::shared_ptr<Sun::Track> Sun::Track::Deserialize(hid_t group){
  double x_   =readDoubleAttribute(group,"x");
  double xini_=readDoubleAttribute(group,"xini");
  double xend_=readDoubleAttribute(group,"xend");
  return std::make_shared<Sun::Track>(x_,xini_,xend_);
}

/*
----------------------------------------------------------------------
         SUN ASNU CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
SunASnu::SunASnu():SunASnu(getResourcePath()+"/astro/bs05_agsop.dat")
{}

SunASnu::SunASnu(std::string sunlocation):Body()
{
  radius = 694439.0*param.km;

  sun_model = quickread(sunlocation);
  arraysize = sun_model.extent(0);

  sun_radius.resize(arraysize);
  sun_density.resize(arraysize);
  sun_xh.resize(arraysize);

  for (unsigned int i=0; i < arraysize;i++){
    sun_radius[i] = sun_model[i][1];
    sun_density[i] = sun_model[i][3];
    sun_xh[i] = sun_model[i][6];
  }

  inter_density=AkimaSpline(sun_radius,sun_density);
  inter_xh=AkimaSpline(sun_radius,sun_xh);
}

SunASnu::SunASnu(std::vector<double> x,std::vector<double> rho,std::vector<double> xh):
Body(),sun_radius(std::move(x)),sun_density(std::move(rho)),sun_xh(std::move(xh)),
inter_density(sun_radius,sun_density),inter_xh(sun_radius,sun_xh)
{
  radius = 695980.0*param.km;
  arraysize = sun_radius.size();
}

// track constructor
SunASnu::Track::Track(double x,double xini,double b_impact):
  Body::Track(x,xini,xini),
  radius_nu(694439.0*param.km),
  b_impact(b_impact)
{
  xend = 2.0*sqrt(SQR(radius_nu)-SQR(b_impact));
}

void SunASnu::Track::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addDoubleAttribute(group,"x",x);
  addDoubleAttribute(group,"xini",xini);
  addDoubleAttribute(group,"xend",xend);
}

std::shared_ptr<SunASnu::Track> SunASnu::Track::Deserialize(hid_t group){
  double x_   =readDoubleAttribute(group,"x");
  double xini_=readDoubleAttribute(group,"xini");
  double xend_=readDoubleAttribute(group,"xend");
  return std::make_shared<SunASnu::Track>(x_,xini_,xend_);
}

void SunASnu::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addUIntAttribute(group,"arraysize",arraysize);
  std::vector<hsize_t> dims {arraysize};
  H5LTmake_dataset_double(group, "sun_radius", 1, dims.data(), sun_radius.data());
  H5LTmake_dataset_double(group, "sun_density", 1, dims.data(), sun_density.data());
  H5LTmake_dataset_double(group, "sun_xh", 1, dims.data(), sun_xh.data());
}

std::shared_ptr<SunASnu> SunASnu::Deserialize(hid_t group){
  unsigned int asize=readUIntAttribute(group,"arraysize");
  std::vector<double> x_vec(asize),rho_vec(asize),xh_vec(asize);
  H5LTread_dataset_double(group,"sun_radius",x_vec.data());
  H5LTread_dataset_double(group,"sun_density",rho_vec.data());
  H5LTread_dataset_double(group,"sun_xh",xh_vec.data());
  return std::make_shared<SunASnu>(x_vec,rho_vec,xh_vec);
}

void SunASnu::Track::FillDerivedParams(std::vector<double>& TrackParams) const{
  TrackParams.push_back(b_impact);
}

double SunASnu::rdensity(double x) const{
// x is adimentional radius : x = 0 : center, x = 1 : radius
  if (x < sun_radius[0]){
    return sun_density[0];
  } else if ( x > sun_radius[arraysize-1]){
    return 0;
  } else {
    return inter_density(x);
  }
}

double SunASnu::rxh(double x) const{
// x is adimentional radius : x = 0 : center, x = 1 : radius
  if (x < sun_radius[0]){
    return sun_xh[0];
  } else if ( x > sun_radius[arraysize-1]){
    return 0;
  } else {
    return inter_xh(x);
  }
}

double SunASnu::density(const GenericTrack& track_input) const
{
  const SunASnu::Track& track_sunasnu = static_cast<const SunASnu::Track&>(track_input);
  double x = track_sunasnu.GetX();
  double b = track_sunasnu.b_impact;

  double r = sqrt(SQR(radius)+SQR(x)-2.0*x*sqrt(SQR(radius)-SQR(b)))/radius;

  return rdensity(r);
}

double SunASnu::ye(const GenericTrack& track_input) const
{
  const SunASnu::Track& track_sunasnu = static_cast<const SunASnu::Track&>(track_input);
  double x = track_sunasnu.GetX();
  double b = track_sunasnu.b_impact;
  double r = sqrt(SQR(radius)+SQR(x)-2.0*x*sqrt(SQR(radius)-SQR(b)))/radius;
  return 0.5*(1.0+rxh(r));
}

SunASnu::~SunASnu(){}

/*
----------------------------------------------------------------------
         EARTHATM CLASS DEFINITIONS
----------------------------------------------------------------------
*/

// constructor
EarthAtm::EarthAtm():EarthAtm(getResourcePath()+"/astro/EARTH_MODEL_PREM_wIso.dat")
{}

// track constructor
EarthAtm::Track::Track(double phi_input,double earth_radius_,double atmheight_):Body::Track(0,0)
{
  earth_radius = earth_radius_;
  atmheight = atmheight_;

  cosphi = cos(phi_input);
  double sinsqphi = 1-cosphi*cosphi;

  L = sqrt(SQR(earth_radius+atmheight)-earth_radius*earth_radius*sinsqphi)-earth_radius*cosphi;

  x = 0.0;
  xini = 0.0;
  xend = L;

  #ifdef EarthAtm_DEBUG
      std::cout << "== Init Track ==" << std::endl;
      std::cout << " phi = " << phi <<
      ", cos(phi) = " << cosphi <<
      ", L = " << radius_nu/param.km << std::endl;
      std::cout << "==" << std::endl;
  #endif
}

void EarthAtm::Track::Serialize(hid_t group) const {
  addH5Attribute(group,"name", GetName());
  addH5Attribute(group,"x",x);
  addH5Attribute(group,"xini",xini);
  addH5Attribute(group,"xend",xend);
  addH5Attribute(group,"cosphi",cosphi);
  addH5Attribute(group,"earth_radius",earth_radius);
  addH5Attribute(group,"atmheight",atmheight);
}

std::shared_ptr<EarthAtm::Track> EarthAtm::Track::Deserialize(hid_t group){
  double x, cosphi;
  readH5Attribute(group,"x", x);
  readH5Attribute(group,"cosphi",cosphi);
  // these parameters were originally not serialized, so give them default values equal to what
  // was originally hard-coded
  double earth_radius = 6371.0*param.km;
  double atmheight = 22.*param.km;
  if(H5Aexists(group,"earth_radius"))
    readH5Attribute(group,"earth_radius", earth_radius);
  if(H5Aexists(group,"atmheight"))
    readH5Attribute(group,"atmheight", atmheight);
  auto track=std::make_shared<EarthAtm::Track>(EarthAtm::Track::MakeWithCosine(cosphi,earth_radius,atmheight));
  track->SetX(x);
  return track;
}

void EarthAtm::Serialize(hid_t group) const {
  addStringAttribute(group,"name", GetName().c_str());
  addUIntAttribute(group,"arraysize",arraysize);
  addUIntAttribute(group,"n_composition",n_composition);
  std::vector<hsize_t> dims {arraysize};
  H5LTmake_dataset_double(group, "earth_radius", 1, dims.data(), earth_radius.data());
  H5LTmake_dataset_double(group, "earth_density", 1, dims.data(), earth_density.data());
  H5LTmake_dataset_double(group, "earth_ye", 1, dims.data(), earth_ye.data());

  // earth_composition is size [n_composition][arraysize] so we must flatten it
  // we will reconstruct it in Deserialize.
  std::vector<double> flattened_earth_composition;
  for (const std::vector<double>& row : earth_composition) {
    flattened_earth_composition.insert(flattened_earth_composition.end(), row.begin(), row.end());
  }
  std::vector<hsize_t> composition_dims = {n_composition, arraysize};
  H5LTmake_dataset_double(group, "earth_composition", 2, composition_dims.data(), flattened_earth_composition.data());

  addH5Attribute(group,"radius",radius);
  addH5Attribute(group,"atm_height",atm_height);
}

std::shared_ptr<EarthAtm> EarthAtm::Deserialize(hid_t group){
  unsigned int asize=readUIntAttribute(group,"arraysize");
  unsigned int ncomposition = 0;
  if(H5Aexists(group,"n_composition"))
    ncomposition = readUIntAttribute(group,"n_composition");
  std::vector<double> x_vec(asize),rho_vec(asize),ye_vec(asize);
  H5LTread_dataset_double(group,"earth_radius",x_vec.data());
  H5LTread_dataset_double(group,"earth_density",rho_vec.data());
  H5LTread_dataset_double(group,"earth_ye",ye_vec.data());
  // these parameters were originally not serialized, so give them default values equal to what
  // was originally hard-coded
  std::vector<std::vector<double>> composition_vec(ncomposition, std::vector<double>(asize));
  if(H5LTfind_dataset(group,"earth_composition")) {
      std::vector<double> flat_composition_vec(ncomposition*asize);
      H5LTread_dataset_double(group,"earth_composition",flat_composition_vec.data());
      for (size_t i = 0; i < ncomposition; ++i) {
          std::copy(flat_composition_vec.begin() + i*asize, flat_composition_vec.begin() + (i+1)*asize, composition_vec[i].begin());
      }
  }

  double radius = 6371.0;
  double atm_height = 22.;
  if(H5Aexists(group,"radius"))
    readH5Attribute(group,"radius", radius);
  if(H5Aexists(group,"atm_height"))
    readH5Attribute(group,"atm_height", atm_height);
  auto earthAtm = std::make_shared<EarthAtm>(x_vec,rho_vec,ye_vec,composition_vec);
  earthAtm->radius = radius;
  earthAtm->SetAtmosphereHeight(atm_height);
  return earthAtm;
}

EarthAtm::Track::Track():
Body::Track(0,0),earth_radius(6371.0*param.km),atmheight(22.*param.km){}

EarthAtm::Track
EarthAtm::Track::MakeWithCosine(double cosphi, double earth_radius_, double atmheight_){
  Track track;

  track.cosphi = cosphi;
  double sinsqphi = 1-track.cosphi*track.cosphi;

  track.L = sqrt(SQR(earth_radius_+atmheight_)-earth_radius_*earth_radius_*sinsqphi)-earth_radius_*cosphi;
  track.x = 0.0;
  track.xini = 0.0;
  track.xend = track.L;

  return(track);
}

void EarthAtm::Track::FillDerivedParams(std::vector<double>& TrackParams) const{
  TrackParams.push_back(acos(cosphi));
  TrackParams.push_back(earth_radius);
  TrackParams.push_back(atmheight);
}
    
EarthAtm::Track EarthAtm::MakeTrack(double phi){
  return Track(phi, GetRadius()*param.km, GetAtmosphereHeight()*param.km);
}
    
EarthAtm::Track EarthAtm::MakeTrackWithCosine(double cosphi){
  return Track::MakeWithCosine(cosphi, GetRadius()*param.km, GetAtmosphereHeight()*param.km);
}
    
double EarthAtm::density(const GenericTrack& track_input) const
{
  const EarthAtm::Track& track_earthatm = static_cast<const EarthAtm::Track&>(track_input);
  double xkm = track_earthatm.GetX()/param.km;
  double sinsqphi = 1-track_earthatm.cosphi*track_earthatm.cosphi;
  double dL = sqrt(SQR(earth_with_atm_radius)-radius*radius*sinsqphi)+radius*track_earthatm.cosphi;
  double r2 = SQR(earth_with_atm_radius) + SQR(xkm) - (track_earthatm.L/param.km+dL)*xkm;
  double r = (r2>0 ? sqrt(r2) : 0);
  
  double rel_r = r/earth_with_atm_radius;
  if ( rel_r < x_radius_min ){
    return x_rho_min;
  }
  else if ( rel_r > x_radius_max and r < radius) {
    return x_rho_max;
  }
  else if ( r > radius) {
    double h = r - radius; // height above surface in km
    double h0 = 7.594; //km obtained by fitting the NRLMSISE atmospheric model up to 60 km to an exponetial profile
    // deviations from NRLMSISE model above this height can be subtantial.
    return 0.0012*exp(-h/h0); // use as constant the atmospheric density at surface in gr/cm^3
  } else {
    return inter_density(rel_r);
  }
}

double EarthAtm::ye(const GenericTrack& track_input) const
{
  const EarthAtm::Track& track_earthatm = static_cast<const EarthAtm::Track&>(track_input);
  double xkm = track_earthatm.GetX()/param.km;
  double sinsqphi = 1-track_earthatm.cosphi*track_earthatm.cosphi;
  double dL = sqrt(SQR(earth_with_atm_radius)-radius*radius*sinsqphi)+radius*track_earthatm.cosphi;
  double r2 = SQR(earth_with_atm_radius) + SQR(xkm) - (track_earthatm.L/param.km+dL)*xkm;
  double r = (r2>0 ? sqrt(r2) : 0);

  double rel_r = r/earth_with_atm_radius;
  if ( rel_r < x_radius_min ){
    return x_ye_min;
  }
  else if ( rel_r > x_radius_max and r < radius) {
    return x_ye_max;
  }
  else if ( r > radius ){
    return 0.494;
  }else {
    return inter_ye(rel_r);
  }
}

std::map<PDGCode, double> EarthAtm::composition(const GenericTrack& track_input) const {
  const EarthAtm::Track& track_earthatm = static_cast<const EarthAtm::Track&>(track_input);
  double xkm = track_earthatm.GetX()/param.km;
  double sinsqphi = 1-track_earthatm.cosphi*track_earthatm.cosphi;
  double dL = sqrt(SQR(earth_with_atm_radius)-radius*radius*sinsqphi)+radius*track_earthatm.cosphi;
  double r2 = SQR(earth_with_atm_radius) + SQR(xkm) - (track_earthatm.L/param.km+dL)*xkm;
  double r = (r2>0 ? sqrt(r2) : 0);

  double rel_r = r/earth_with_atm_radius;
  if ( rel_r < x_radius_min ){
    return x_composition_min;
  }
  else if ( rel_r > x_radius_max and r < radius ) {
    return x_composition_max;
  }
  else if ( r > radius ){
    return x_composition_max; // TODO: Check if this is right -PW
  } else {
    std::map<PDGCode, double> composition_map;
    for (const auto& pair : inter_composition) {
        composition_map.emplace(pair.first, pair.second(rel_r));
    }
    return composition_map;
  }
}

EarthAtm::EarthAtm(std::string filepath):Body()
{
  // Note this radius is the averaged Earth radius. -PW
  // Equitorial radius = 6378.1370 km
  // Polar radius = 6356.7523 km
  radius = 6371.0; // km
  atm_height = 22; // km
  earth_with_atm_radius = radius + atm_height;

  marray<double,2> earth_model = quickread(filepath);
  arraysize = earth_model.extent(0);

  // Assume the first 3 columns of the file are the regular PREM details:
  // r/R, rho, ye, [nuclear target fractions]
  // Maybe we should check that a "good" file is being used, i.e. has the required columns -PW
  n_composition = earth_model.extent(1) - 3;

  earth_radius.resize(arraysize);
  earth_density.resize(arraysize);
  earth_ye.resize(arraysize);
  earth_composition.resize(n_composition); // fraction of each component
  for (std::vector<double>& composition_vector : earth_composition) {
    composition_vector.resize(arraysize);
  }

  for (size_t i = 0; i < arraysize; i++){
    earth_radius[i] = earth_model[i][0];
    earth_density[i] = earth_model[i][1];
    earth_ye[i] = earth_model[i][2];
    for (size_t j = 0; j < n_composition; j++) {
      earth_composition[j][i] = earth_model[i][3+j];
    }
  }

  x_radius_min = earth_radius[0];
  x_radius_max = earth_radius[arraysize-1];
  x_rho_min = earth_density[0];
  x_rho_max = earth_density[arraysize-1];
  x_ye_min = earth_ye[0];
  x_ye_max = earth_ye[arraysize-1];

  // The PREM data file has no meta data, so we must assume you are using these elements
  // TODO: consider a simpler model that only consists of a few elements -PW
  std::vector<PDGCode> composition_codes = { proton, oxygen, sodium, magnesium, aluminum, silicon, sulfur, calcium, iron, nickel };
	
  inter_density=AkimaSpline(earth_radius,earth_density);
  inter_ye=AkimaSpline(earth_radius,earth_ye);
  for (size_t i = 0; i < n_composition; i++) {
    PDGCode tgt_id = composition_codes[i];
    inter_composition[tgt_id] = AkimaSpline(earth_radius, earth_composition[i]);
    x_composition_min[tgt_id] = earth_composition[i][0];
    x_composition_max[tgt_id] = earth_composition[i][arraysize-1];
  }
}

EarthAtm::EarthAtm(std::vector<double> x,std::vector<double> rho,std::vector<double> ye, std::vector<std::vector<double>> composition):
Body(),earth_radius(std::move(x)),earth_density(std::move(rho)),earth_ye(std::move(ye)),
earth_composition(composition),inter_density(earth_radius,earth_density),inter_ye(earth_radius,earth_ye)
{
  assert("nuSQUIDS::Error::EarthConstructor: Invalid array sizes." && earth_radius.size() == earth_density.size() && earth_radius.size() == earth_ye.size());
  // The Input file should have the radius specified from 0 to 1.
  // where 0 is the center of the Earth and 1 is the surface.
  radius = 6371.0; // km
  atm_height = 22; // km
  earth_with_atm_radius = radius + atm_height;
  arraysize = earth_radius.size();
  n_composition = earth_composition.size();

  x_radius_min = earth_radius[0];
  x_radius_max = earth_radius[arraysize-1];
  x_rho_min = earth_density[0];
  x_rho_max = earth_density[arraysize-1];
  x_ye_min = earth_ye[0];
  x_ye_max = earth_ye[arraysize-1];

  // TODO: serialize this?
  std::vector<PDGCode> composition_codes = { proton, oxygen, sodium, magnesium, aluminum, silicon, sulfur, calcium, iron, nickel };
  for (size_t i = 0; i < n_composition; i++) {
    PDGCode tgt_id = composition_codes[i];
    inter_composition[tgt_id] = AkimaSpline(earth_radius, earth_composition[i]);
    x_composition_min[tgt_id] = earth_composition[i][0];
    x_composition_max[tgt_id] = earth_composition[i][arraysize-1];
  }
}

EarthAtm::~EarthAtm(){}
    
void EarthAtm::SetAtmosphereHeight(double height){
  atm_height = height;
  earth_with_atm_radius = radius + atm_height;
}

// body registration stuff

std::map<std::string,std::function<std::shared_ptr<Body>(hid_t)>>* body_registry=NULL;
std::map<std::string,std::function<std::shared_ptr<Track>(hid_t)>>* track_registry=NULL;

namespace detail{
  registerBody::registerBody(const std::string& name,std::function<std::shared_ptr<Body>(hid_t)> fdeserialize){
    if(!body_registry)
      body_registry=new std::map<std::string,std::function<std::shared_ptr<Body>(hid_t)>>;
    body_registry->insert(std::make_pair(name,fdeserialize));
  }
  registerTrack::registerTrack(const std::string& name,std::function<std::shared_ptr<Track>(hid_t)> fdeserialize){
    if(!track_registry)
      track_registry=new std::map<std::string,std::function<std::shared_ptr<Track>(hid_t)>>;
    track_registry->insert(std::make_pair(name,fdeserialize));
  }
}

std::function<std::shared_ptr<Body>(hid_t)> GetBodyDeserializer(std::string body_name){
  auto it=body_registry->find(body_name);
  if(it==body_registry->end())
    throw std::runtime_error("Unknown Body type: "+body_name);
  return it->second;
}

std::function<std::shared_ptr<Track>(hid_t)> GetTrackDeserializer(std::string track_name){
  auto it=track_registry->find(track_name);
  if(it==track_registry->end())
    throw std::runtime_error("Unknown track type: "+track_name);
  return it->second;
}

} // close namespace

// registering the default bodies
using namespace nusquids;
NUSQUIDS_REGISTER_BODY(Vacuum);
NUSQUIDS_REGISTER_BODY(ConstantDensity);
NUSQUIDS_REGISTER_BODY(VariableDensity);
NUSQUIDS_REGISTER_BODY(Earth);
NUSQUIDS_REGISTER_BODY(EarthAtm);
NUSQUIDS_REGISTER_BODY(Sun);
NUSQUIDS_REGISTER_BODY(SunASnu);
