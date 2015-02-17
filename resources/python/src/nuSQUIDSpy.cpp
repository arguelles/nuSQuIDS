#include <boost/python.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/to_python_converter.hpp>
#include "container_conversions.h"
#include <SQuIDS/SQUIDS.h>
#include <nuSQuIDS/nuSQUIDS.h>
#include <nuSQuIDS/marray.h>

#include <numpy/ndarrayobject.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>


using namespace boost::python;
using namespace nusquids;
namespace bp = boost::python;

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


// converting marray to numpy array and back
template<unsigned int DIM>
static boost::python::object marray_to_numpyarray( marray<double,DIM> const & iarray){
  // get the data from the marray
  double * data = iarray.size() ? const_cast<double*>(iarray.get_data()) : static_cast<double*>(NULL);
  // construct numpy object
  npy_intp * size = new npy_intp[DIM];
  for(unsigned int i = 0; i < DIM; i++)
    size[i] = iarray.extent(i);
  PyObject * pyObj = PyArray_SimpleNewFromData(DIM,size,NPY_DOUBLE,data);
  boost::python::handle<> handle(pyObj);
  boost::python::numeric::array arr(handle);

  delete [] size;

  // return numpy object
  return arr.copy();
}

template<typename T,unsigned int DIM>
static marray<T,DIM> numpyarray_to_marray(PyObject * iarray, NPY_TYPES type_num){
  // es un array de numpy
  if (! PyArray_Check(iarray) )
  {
    PyErr_SetString(PyExc_TypeError, "numpyarray_to_marray: Input is not a numpy array.");
    boost::python::throw_error_already_set();
  }
  // si es que fuera un array de numpy castearlo
  //PyArrayObject* numpy_array = (PyArrayObject*) iarray;
  // lets get the contiguos C-style array
  PyArrayObject* numpy_array = PyArray_GETCONTIGUOUS((PyArrayObject*)iarray);

  // revisemos que los tipos del array sean dobles o que
  if ( PyArray_DESCR(numpy_array)->type_num != type_num )
  {
    if ( PyArray_DESCR(numpy_array)->type_num == NPY_LONG &&
        PyArray_ITEMSIZE(numpy_array) == 4 && type_num == NPY_INT)
    {
      // numpy on 32 bits sets numpy.int32 to NPY_LONG. So its all ok.
    }
    else
    {
      PyErr_SetString(PyExc_TypeError, "numpyarray_to_marray: numpy type is not the same as the input array type.");
      boost::python::throw_error_already_set();
    }
  }

  // arrays vacios
  if (PyArray_SIZE(numpy_array) == 0){
      PyErr_SetString(PyExc_TypeError,"numpyarray_to_marray: empty numpy array.");
      boost::python::throw_error_already_set();
  }

  // create numpy iterator
  NpyIter* iter = NpyIter_New(numpy_array, NPY_ITER_READONLY|
                             NPY_ITER_EXTERNAL_LOOP|
                             NPY_ITER_REFS_OK,
                             NPY_KEEPORDER, NPY_NO_CASTING,
                             NULL);

  unsigned int array_dim = PyArray_NDIM(numpy_array);
  assert(DIM == array_dim && "No matching dimensions.");

  // get numpy array shape and create marray object
  npy_intp* array_shape = PyArray_SHAPE(numpy_array);
  std::vector<size_t> dimensions;
  for(unsigned int i = 0; i < array_dim; i++)
    dimensions.push_back(array_shape[i]);

  // construct output object
  marray<T,DIM> oarray;
  oarray.resize(dimensions);
  auto it = oarray.begin();

  NpyIter_IterNextFunc *iternext = NpyIter_GetIterNext(iter, NULL);
  char** dataptr = NpyIter_GetDataPtrArray(iter);
  npy_intp* strideptr = NpyIter_GetInnerStrideArray(iter);
  npy_intp* sizeptr = NpyIter_GetInnerLoopSizePtr(iter);
  npy_intp iop, nop = NpyIter_GetNOp(iter);

  do{
    char* data = *dataptr;
    npy_intp count = *sizeptr;
    npy_intp stride = *strideptr;

    while (count--)
    {
      for (iop = 0; iop < nop; ++iop, data+=stride)
        *it++ = *(T*)data;
    }
  } while(iternext(iter));

  NpyIter_Deallocate(iter);

  return oarray;
}

// nuSQUIDS wrap functions
static void wrap_WriteStateHDF5(nuSQUIDS* nusq, std::string path){
  nusq->WriteStateHDF5(path);
}

static void wrap_ReadStateHDF5(nuSQUIDS* nusq, std::string path){
  nusq->ReadStateHDF5(path);
}

static void wrap_Set_initial_state(nuSQUIDS* nusq, PyObject * array, std::string neutype){
  if (! PyArray_Check(array) )
  {
    throw std::runtime_error("nuSQUIDSpy::Error:Input array is not a numpy array.");
  }
  PyArrayObject* numpy_array = (PyArrayObject*)array;
  unsigned int array_dim = PyArray_NDIM(numpy_array);


  if ( array_dim == 1 ) {
    marray<double,1> state = numpyarray_to_marray<double,1>(array, NPY_DOUBLE);
    nusq->Set_initial_state(state,neutype);
  } else if ( array_dim == 2 ) {
    marray<double,2> state = numpyarray_to_marray<double,2>(array, NPY_DOUBLE);
    nusq->Set_initial_state(state,neutype);
  } else if ( array_dim == 3 ) {
    marray<double,3> state = numpyarray_to_marray<double,3>(array, NPY_DOUBLE);
    nusq->Set_initial_state(state,neutype);
  } else
    throw std::runtime_error("nuSQUIDS::Error:Input array has wrong dimenions.");
}

static void wrap_Set_initial_state_atm(nuSQUIDSAtm* nusq_atm, PyObject * array, std::string neutype){
  /*
  bool isint = false;
  //if ( array.get_dtype() == np::dtype::get_builtin<int>() | )
  // int64
  //  isint = true;
  if ( array.get_dtype() != np::dtype::get_builtin<double>() )
    isint = true;
    //throw std::runtime_error("nuSQUIDS::Input array cannot be converted to double.");

  Py_intptr_t const * strides = array.get_strides();


  if ( array.get_nd() == 3 ) {
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
    nusq_atm->Set_initial_state(state,neutype);
   } else if ( array.get_nd() == 4 ) {
    std::vector < std::vector< std::vector < std::vector<double> > > >  state(array.shape(0));
    for (int i = 0; i < array.shape(0); i++){
      state[i].resize(array.shape(1));
      for (int j = 0; j < array.shape(1); j++){
        state[i][j].resize(array.shape(2));
        for (int k = 0; k < array.shape(2); k++){
          state[i][j][k].resize(array.shape(3));
          for (int l = 0; l < array.shape(3); l++){
            if (isint) {
              state[i][j][k][l] = (double)*reinterpret_cast<const int*>(array.get_data() + i*strides[0] + j*strides[1] + k*strides[2] + l*strides[3]);
            } else {
              state[i][j][k][l] = (double)*reinterpret_cast<const double *>(array.get_data() + i*strides[0] + j*strides[1] + k*strides[2] + l*strides[3]);
              //std::cout << i << " " << j << " " << k << " " << state[i][j][k] << std::endl;
            }
          }
        }
      }
    }
    nusq_atm->Set_initial_state(state,neutype);
  } else
    throw std::runtime_error("nuSQUIDSAtm::Error:Input array has wrong dimenions.");
*/
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

  // import numpy array definitions
  import_array();
  import_ufunc();

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
    .def(init<unsigned int>())
    .def("Rotate",&SU_vector::Rotate)
    .def("Dim",&SU_vector::Dim)
    .def("GetComponents",&SU_vector::GetComponents)
  ;

  enum_<NeutrinoType>("NeutrinoType")
    .value("neutrino",neutrino)
    .value("antineutrino",antineutrino)
    .value("both",both)
  ;

  class_<nuSQUIDS, boost::noncopyable, std::shared_ptr<nuSQUIDS> >("nuSQUIDS", init<double,double,unsigned int,unsigned int,NeutrinoType,bool,bool>())
    .def(init<std::string>())
    .def(init<unsigned int,NeutrinoType>())
    .def("Set_initial_state",wrap_Set_initial_state)
    .def("Set_Body",&nuSQUIDS::Set_Body, bp::arg("Body"))
    .def("Set_Track",&nuSQUIDS::Set_Track, bp::arg("Track"))
    .def("Set_E",&nuSQUIDS::Set_E, bp::arg("NeutrinoEnergy"))
    .def("EvolveState",&nuSQUIDS::EvolveState)
    .def("GetERange",&nuSQUIDS::GetERange)
    .def("WriteStateHDF5",&nuSQUIDS::WriteStateHDF5)
    .def("WriteStateHDF5",wrap_WriteStateHDF5)
    .def("ReadStateHDF5",&nuSQUIDS::ReadStateHDF5)
    .def("ReadStateHDF5",wrap_ReadStateHDF5)
    .def("H0",&nuSQUIDS::H0)
    .def("HI",(SU_vector(nuSQUIDS::*)(unsigned int,unsigned int,double) const)&nuSQUIDS::HI)
    .def("HI",(SU_vector(nuSQUIDS::*)(unsigned int,unsigned int) const)&nuSQUIDS::HI)
    .def("GetNumNeu",&nuSQUIDS::GetNumNeu)
    .def("EvalMass",(double(nuSQUIDS::*)(unsigned int) const)&nuSQUIDS::EvalMass)
    .def("EvalFlavor",(double(nuSQUIDS::*)(unsigned int) const)&nuSQUIDS::EvalFlavor)
    .def("EvalMass",(double(nuSQUIDS::*)(unsigned int,double,unsigned int) const)&nuSQUIDS::EvalMass)
    .def("EvalFlavor",(double(nuSQUIDS::*)(unsigned int,double,unsigned int) const)&nuSQUIDS::EvalFlavor)
    .def("EvalMassAtNode",(double(nuSQUIDS::*)(unsigned int,unsigned int,unsigned int) const)&nuSQUIDS::EvalMassAtNode)
    .def("EvalFlavorAtNode",(double(nuSQUIDS::*)(unsigned int,unsigned int,unsigned int) const)&nuSQUIDS::EvalFlavorAtNode)
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

  class_<nuSQUIDSAtm, boost::noncopyable, std::shared_ptr<nuSQUIDSAtm> >("nuSQUIDSAtm", init<double,double,unsigned int,double,double,unsigned int,unsigned int,NeutrinoType,bool,bool>())
    .def(init<std::string>())
    .def("EvolveState",&nuSQUIDSAtm::EvolveState)
    .def("Set_TauRegeneration",&nuSQUIDSAtm::Set_TauRegeneration)
    .def("EvalFlavor",&nuSQUIDSAtm::EvalFlavor)
    .def("WriteStateHDF5",&nuSQUIDSAtm::WriteStateHDF5)
    .def("ReadStateHDF5",&nuSQUIDSAtm::ReadStateHDF5)
    .def("Set",&nuSQUIDSAtm::Set)
    .def("Set_ProgressBar",&nuSQUIDSAtm::Set_ProgressBar)
    .def("Set_MixingParametersToDefault",&nuSQUIDSAtm::Set_MixingParametersToDefault)
    .def_readonly("units", &nuSQUIDSAtm::units)
    .def("Set_rel_error",&nuSQUIDSAtm::Set_rel_error)
    .def("Set_abs_error",&nuSQUIDSAtm::Set_abs_error)
    .def("GetNumE",&nuSQUIDSAtm::GetNumE)
    .def("GetNumCos",&nuSQUIDSAtm::GetNumCos)
    .def("Set_initial_state",wrap_Set_initial_state_atm)
    .def("GetERange",&nuSQUIDSAtm::GetERange)
    .def("GetCosthRange",&nuSQUIDSAtm::GetCosthRange)
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
    .def(init<std::string>())
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
    .def(init<std::string>())
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
