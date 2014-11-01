#include <boost/python.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/numpy.hpp>
#include "container_conversions.h"
#include "nuSQUIDS.h"
#include "SQUIDS.h"

using namespace boost::python;
using namespace nusquids;
namespace bp = boost::python;
namespace np = boost::numpy;

template<class T>
struct VecToList
{
  static PyObject* convert(const std::vector<T>& vec){
    boost::python::list* l = new boost::python::list();
    for(size_t i =0; i < vec.size(); i++)
      (*l).append(vec[i]);

    return l->ptr();
  }
};

// nuSQUIDS wrap functions

static void wrap_Set_initial_state(nuSQUIDS* nusq, const np::ndarray & array, std::string neutype){
  bool isint = false;
  //if ( array.get_dtype() == np::dtype::get_builtin<int>() | )
  // int64
  //  isint = true;
  if ( array.get_dtype() != np::dtype::get_builtin<double>() )
    isint = true;
    //throw std::runtime_error("nuSQUIDS::Input array cannot be converted to double.");

  Py_intptr_t const * strides = array.get_strides();
  if ( array.get_nd() == 1 ) {
    std::vector<double> state(array.shape(0));
    for (int i = 0; i < array.shape(0); i++){
      if (isint)
        state[i] = (double)*reinterpret_cast<const int*>(array.get_data() + i*strides[0]);
      else
        state[i] = (double)*reinterpret_cast<const double *>(array.get_data() + i*strides[0]);
    }
    nusq->Set_initial_state(state,neutype);
  } else if ( array.get_nd() == 2 ) {
    std::vector< std::vector<double> > state(array.shape(0));
    for (int i = 0; i < array.shape(0); i++){
      state[i].resize(array.shape(1));
      for (int j = 0; j < array.shape(1); j++){
        if (isint)
          state[i][j] = (double)*reinterpret_cast<const int*>(array.get_data() + i*strides[0] + j*strides[1]);
        else
          state[i][j] = (double)*reinterpret_cast<const double *>(array.get_data() + i*strides[0] + j*strides[1]);
      }
    }
    nusq->Set_initial_state(state,neutype);
  } else if ( array.get_nd() == 3 ) {
    std::vector< std::vector < std::vector<double> > > state(array.shape(0));
    for (int i = 0; i < array.shape(0); i++){
      state[i].resize(array.shape(1));
      for (int j = 0; j < array.shape(1); j++){
        state[i][j].resize(array.shape(2));
        for (int k = 0; k < array.shape(2); k++){
          if (isint) {
            state[i][j][k] = (double)*reinterpret_cast<const int*>(array.get_data() + i*strides[0] + j*strides[1] + k*strides[2]);
          } else {
            state[i][j][k] = (double)*reinterpret_cast<const double *>(array.get_data() + i*strides[0] + j*strides[1] + k*strides[2]);
            //std::cout << i << " " << j << " " << k << " " << state[i][j][k] << std::endl;
          }
        }
      }
    }
    nusq->Set_initial_state(state,neutype);
  } else
    throw std::runtime_error("nuSQUIDS::Error:Input array has wrong dimenions.");
}

enum GSL_STEP_FUNCTIONS {
  GSL_STEP_RK2,
  GSL_STEP_RK4,
  GSL_STEP_RKF45,
  GSL_STEP_RKCK,
  GSL_STEP_RK8PD,
  /*
  GSL_STEP_RK1IMP,
  GSL_STEP_RK2IMP,
  GSL_STEP_RK4IMP,
  GSL_STEP_BSIMP,
  GSL_STEP_MSBDF,
  */
  GSL_STEP_MSADAMS
};

static void wrap_Set_GSL_STEP(nuSQUIDS* nusq, GSL_STEP_FUNCTIONS step_enum){
  switch(step_enum){
    case GSL_STEP_RK2:
      nusq->Set_GSL_step(gsl_odeiv2_step_rk2);
      break;
    case GSL_STEP_RK4:
      nusq->Set_GSL_step(gsl_odeiv2_step_rk4);
      break;
    case GSL_STEP_RKF45:
      nusq->Set_GSL_step(gsl_odeiv2_step_rkf45);
      break;
    case GSL_STEP_RKCK:
      nusq->Set_GSL_step(gsl_odeiv2_step_rkck);
      break;
    case GSL_STEP_RK8PD:
      nusq->Set_GSL_step(gsl_odeiv2_step_rk8pd);
      break;
      /*
    case GSL_STEP_RK1IMP:
      nusq->Set_GSL_step(gsl_odeiv2_step_rk1imp);
      break;
    case GSL_STEP_RK2IMP:
      nusq->Set_GSL_step(gsl_odeiv2_step_rk2imp);
      break;
    case GSL_STEP_RK4IMP:
      nusq->Set_GSL_step(gsl_odeiv2_step_rk4imp);
      break;
    case GSL_STEP_BSIMP:
      nusq->Set_GSL_step(gsl_odeiv2_step_bsimp);
      break;
    case GSL_STEP_MSBDF:
      nusq->Set_GSL_step(gsl_odeiv2_step_msbdf);
      break;
      */
    case GSL_STEP_MSADAMS:
      nusq->Set_GSL_step(gsl_odeiv2_step_msadams);
      break;
  }
}

// nuSQUIDSpy module definitions

BOOST_PYTHON_MODULE(nuSQUIDSpy)
{
  //Py_Initialize();
  np::initialize();

  enum_<GSL_STEP_FUNCTIONS>("GSL_STEP_FUNCTIONS")
    .value("GSL_STEP_RK2",GSL_STEP_RK2)
    .value("GSL_STEP_RK4",GSL_STEP_RK4)
    .value("GSL_STEP_RKF45",GSL_STEP_RKF45)
    .value("GSL_STEP_RKCK",GSL_STEP_RKCK)
    .value("GSL_STEP_RK8PD",GSL_STEP_RK8PD)
    /*
    .value("GSL_STEP_RK1IMP",GSL_STEP_RK1IMP)
    .value("GSL_STEP_RK2IMP",GSL_STEP_RK2IMP)
    .value("GSL_STEP_RK4IMP",GSL_STEP_RK4IMP)
    .value("GSL_STEP_BSIMP",GSL_STEP_BSIMP)
    .value("GSL_STEP_MSBDF",GSL_STEP_MSBDF)
    */
    .value("GSL_STEP_MSADAMS",GSL_STEP_MSADAMS)
    ;

  enum_<BASIS>("BASIS")
    .value("MASS",mass)
    .value("INTERACTION",interaction)
    ;

  enum_<MixingParameter>("MixingParameter")
    .value("TH12",TH12)
    .value("TH13",TH13)
    .value("TH23",TH23)
    .value("TH14",TH14)
    .value("TH24",TH24)
    .value("TH34",TH34)
    .value("TH15",TH15)
    .value("TH25",TH25)
    .value("TH35",TH35)
    .value("TH45",TH45)
    .value("TH16",TH16)
    .value("TH26",TH26)
    .value("TH36",TH36)
    .value("TH46",TH46)
    .value("TH56",TH56)
    .value("DELTA1",DELTA1)
    .value("DELTA2",DELTA2)
    .value("DELTA3",DELTA3)
    .value("DM21SQ",DM21SQ)
    .value("DM31SQ",DM31SQ)
    .value("DM41SQ",DM41SQ)
    .value("DM51SQ",DM51SQ)
    .value("DM61SQ",DM61SQ)
    ;

  class_<SU_vector, boost::noncopyable,std::shared_ptr<SU_vector> >("SU_vector")
    .def(init< std::vector<double> >())
    .def(init<int>())
    .def("Rescale",&SU_vector::Rescale)
    .def("Rotate",&SU_vector::Rotate)
    .def("Dim",&SU_vector::Dim)
    .def("GetComponents",&SU_vector::GetComponents)
  ;

  //class_<SQUIDS>("SQUIDS")
  //  .def("Set",&SQUIDS::Set)
  //;

  class_<nuSQUIDS, boost::noncopyable, std::shared_ptr<nuSQUIDS> >("nuSQUIDS", init<double,double,int,int,std::string,bool,bool>())
    .def(init<std::string>())
    .def(init<int,std::string>())
    .def("Set_initial_state",
        (void(nuSQUIDS::*)(std::vector<double>,std::string))&nuSQUIDS::Set_initial_state,
        ( bp::arg("InitialState"), bp::arg("NeutrinoType") )
        )
    .def("Set_initial_state",wrap_Set_initial_state)
    .def("Set_Body",&nuSQUIDS::Set_Body, bp::arg("Body"))
    .def("Set_Track",&nuSQUIDS::Set_Track, bp::arg("Track"))
    .def("Set_E",&nuSQUIDS::Set_E, bp::arg("NeutrinoEnergy"))
    .def("EvolveState",&nuSQUIDS::EvolveState)
    .def("GetERange",&nuSQUIDS::GetERange)
    .def("WriteStateHDF5",&nuSQUIDS::WriteStateHDF5)
    .def("ReadStateHDF5",&nuSQUIDS::ReadStateHDF5)
    .def("H0",&nuSQUIDS::H0)
    .def("HI",(SU_vector(nuSQUIDS::*)(int,double))&nuSQUIDS::HI)
    .def("HI",(SU_vector(nuSQUIDS::*)(int))&nuSQUIDS::HI)
    .def("GetNumNeu",&nuSQUIDS::GetNumNeu)
    .def("EvalMass",(double(nuSQUIDS::*)(int))&nuSQUIDS::EvalMass)
    .def("EvalFlavor",(double(nuSQUIDS::*)(int))&nuSQUIDS::EvalFlavor)
    .def("EvalMass",(double(nuSQUIDS::*)(int,double,int))&nuSQUIDS::EvalMass)
    .def("EvalFlavor",(double(nuSQUIDS::*)(int,double,int))&nuSQUIDS::EvalFlavor)
    .def("EvalMassAtNode",(double(nuSQUIDS::*)(int,int,int))&nuSQUIDS::EvalMassAtNode)
    .def("EvalFlavorAtNode",(double(nuSQUIDS::*)(int,int,int))&nuSQUIDS::EvalFlavorAtNode)
    .def("GetERange",&nuSQUIDS::GetERange)
    .def("GetHamiltonian",&nuSQUIDS::GetHamiltonian)
    .def("GetState",&nuSQUIDS::GetState)
//    .def("Set",(void(nuSQUIDS::*)(string,double))&nuSQUIDS::Set)
//    .def("Set",(void(nuSQUIDS::*)(string,bool))&nuSQUIDS::Set)
//    .def("Set",(void(nuSQUIDS::*)(string,int))&nuSQUIDS::Set)
    .def("Set_h_min",&nuSQUIDS::Set_h_min)
    .def("Set_h_max",&nuSQUIDS::Set_h_max)
    .def("Set_h",&nuSQUIDS::Set_h)
    .def("Set_rel_error",&nuSQUIDS::Set_rel_error)
    .def("Set_abs_error",&nuSQUIDS::Set_abs_error)
    .def("Set_AdaptiveStep",&nuSQUIDS::Set_AdaptiveStep)
    .def("Set_GSL_step",wrap_Set_GSL_STEP)
    .def("Set_TauRegeneration",&nuSQUIDS::Set_TauRegeneration)
    .def("Set_ProgressBar",&nuSQUIDS::Set_ProgressBar)
    .def("Set_MixingParametersToDefault",&nuSQUIDS::Set_MixingParametersToDefault)
    .def("Set_Basis",&nuSQUIDS::Set_Basis)
    .def("Set",&nuSQUIDS::Set)
    .def("GetTrack",&nuSQUIDS::GetTrack)
    .def("GetBody",&nuSQUIDS::GetBody)
    .def("GetNumE",&nuSQUIDS::GetNumE)
    .def_readonly("units", &nuSQUIDS::units)
  ;


  class_<Const>("Const")
    .def_readonly("TeV",&Const::TeV)
    .def_readonly("GeV",&Const::GeV)
    .def_readonly("MeV",&Const::MeV)
    .def_readonly("keV",&Const::keV)
    .def_readonly("eV",&Const::eV)
    .def_readonly("kg",&Const::kg)
    .def_readonly("gr",&Const::gr)
    .def_readonly("meter",&Const::meter)
    .def_readonly("cm",&Const::cm)
    .def_readonly("km",&Const::km)
    .def_readonly("fermi",&Const::fermi)
    .def_readonly("angstrom",&Const::angstrom)
    .def_readonly("AU",&Const::AU)
    .def_readonly("parsec",&Const::parsec)
    .def_readonly("pb",&Const::picobarn)
    .def_readonly("fb",&Const::femtobarn)
    .def_readonly("sec",&Const::sec)
    .def_readonly("hour",&Const::hour)
    .def_readonly("day",&Const::day)
    .def_readonly("year",&Const::year)
    .def_readonly("_s",&Const::s)
    .def_readonly("_c",&Const::c)
    .def_readonly("_dcp",&Const::dcp)
  ;

  {
    scope outer
    = class_<Body, std::shared_ptr<Body> >("Body")
    .def("density",&Body::density)
    .def("ye",&Body::ye)
    ;

    class_<Body::Track, std::shared_ptr<Body::Track> >("Track")
    .def("GetInitialX",&Body::Track::GetInitialX)
    .def("GetFinalX",&Body::Track::GetFinalX)
    .def("GetX",&Body::Track::GetX)
    .def("SetX",&Body::Track::SetX)
    ;
  }

  {
    scope outer
    = class_<Vacuum, bases<Body>, std::shared_ptr<Vacuum> >("Vacuum")
    ;

    class_<Vacuum::Track, std::shared_ptr<Vacuum::Track> >("Track")
    .def(init<double,double>())
    .def("density",&Vacuum::density)
    .def("ye",&Vacuum::ye)
    .def("GetInitialX",&Vacuum::Track::GetInitialX)
    .def("GetFinalX",&Vacuum::Track::GetFinalX)
    .def("GetX",&Vacuum::Track::GetX)
    .def("SetX",&Vacuum::Track::SetX)
    ;

    implicitly_convertible< std::shared_ptr<Vacuum>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<Vacuum::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<ConstantDensity, bases<Body>, std::shared_ptr<ConstantDensity> >("ConstantDensity")
    .def(init<double,double>())
    ;

    class_<ConstantDensity::Track, std::shared_ptr<ConstantDensity::Track> >("Track")
    .def(init<double,double>())
    .def("GetInitialX",&ConstantDensity::Track::GetInitialX)
    .def("GetFinalX",&ConstantDensity::Track::GetFinalX)
    .def("GetX",&ConstantDensity::Track::GetX)
    .def("SetX",&ConstantDensity::Track::SetX)
    ;

    implicitly_convertible< std::shared_ptr<ConstantDensity>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<ConstantDensity::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<VariableDensity, bases<Body>, std::shared_ptr<VariableDensity> >("VariableDensity")
    .def(init< std::vector<double>,std::vector<double>,std::vector<double> >())
    ;

    class_<VariableDensity::Track, std::shared_ptr<VariableDensity::Track> >("Track")
    .def(init<double,double>())
    .def("GetInitialX",&VariableDensity::Track::GetInitialX)
    .def("GetFinalX",&VariableDensity::Track::GetFinalX)
    .def("GetX",&VariableDensity::Track::GetX)
    .def("SetX",&VariableDensity::Track::SetX)
    ;

    implicitly_convertible< std::shared_ptr<VariableDensity>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<VariableDensity::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<Earth, bases<Body>, std::shared_ptr<Earth> >("Earth")
    .def(init<string>())
    ;

    class_<Earth::Track, std::shared_ptr<Earth::Track> >("Track")
    .def(init<double,double,double>())
    .def("GetInitialX",&Earth::Track::GetInitialX)
    .def("GetFinalX",&Earth::Track::GetFinalX)
    .def("GetX",&Earth::Track::GetX)
    .def("SetX",&Earth::Track::SetX)
    ;

    implicitly_convertible< std::shared_ptr<Earth>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<Earth::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<Sun, bases<Body>, std::shared_ptr<Sun> >("Sun")
    ;

    class_<Sun::Track, std::shared_ptr<Sun::Track> >("Track")
    .def(init<double,double>())
    .def("GetInitialX",&Sun::Track::GetInitialX)
    .def("GetFinalX",&Sun::Track::GetFinalX)
    .def("GetX",&Sun::Track::GetX)
    .def("SetX",&Sun::Track::SetX)
    ;

    implicitly_convertible< std::shared_ptr<Sun>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<Sun::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<SunASnu, bases<Body>, std::shared_ptr<SunASnu> >("SunASnu")
    ;

    class_<SunASnu::Track, std::shared_ptr<SunASnu::Track> >("Track")
    .def(init<double,double>())
    .def("GetInitialX",&SunASnu::Track::GetInitialX)
    .def("GetFinalX",&SunASnu::Track::GetFinalX)
    .def("GetX",&SunASnu::Track::GetX)
    .def("SetX",&SunASnu::Track::SetX)
    ;

    implicitly_convertible< std::shared_ptr<SunASnu>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<SunASnu::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<EarthAtm, bases<Body>, std::shared_ptr<EarthAtm> >("EarthAtm")
    .def(init<string>())
    ;

    class_<EarthAtm::Track, std::shared_ptr<EarthAtm::Track> >("Track")
    .def(init<double>())
    .def("GetInitialX",&EarthAtm::Track::GetInitialX)
    .def("GetFinalX",&EarthAtm::Track::GetFinalX)
    .def("GetX",&EarthAtm::Track::GetX)
    .def("SetX",&EarthAtm::Track::SetX)
    ;

    implicitly_convertible< std::shared_ptr<EarthAtm>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<EarthAtm::Track>, std::shared_ptr<Body::Track> >();
  }

  // python cotainer to vector<double> convertion
  using namespace scitbx::boost_python::container_conversions;
  from_python_sequence< std::vector<double>, variable_capacity_policy >();
  //from_python_sequence< std::vector< std::vector<double> >, variable_capacity_policy >();
  to_python_converter< std::vector<double, class std::allocator<double> >, VecToList<double> > ();

}
