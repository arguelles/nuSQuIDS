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

#define H5Gopen_vers 2
#define H5Gcreate_vers 2
#define H5Eset_auto_vers 2
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <boost/python.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python/overloads.hpp>
#include "container_conversions.h"
#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>
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
struct marray_to_numpyarray {
  static PyObject* convert( marray<double,DIM> const & iarray){
    // get the data from the marray
    double * data = iarray.size() ? const_cast<double*>(iarray.get_data()) : static_cast<double*>(NULL);
    // construct numpy object
    npy_intp size[DIM];
    for(unsigned int i = 0; i < DIM; i++)
      size[i] = iarray.extent(i);
    PyArrayObject * pyObj = (PyArrayObject*) PyArray_SimpleNew(DIM,size,PyArray_DOUBLE);
    memcpy(pyObj->data, data, sizeof(double) * iarray.size());

    return PyArray_Return(pyObj);
  }
};

template<unsigned int Dim>
struct marray_from_python{
  marray_from_python(){
    boost::python::converter::registry::push_back(&convertible,
                                                  &construct,
                                                  boost::python::type_id<marray<double,Dim>>());
  }
  
  static void* convertible(PyObject* obj_ptr){
    //accept only numpy arrays
    if(!PyArray_Check(obj_ptr))
      return(NULL);
    PyArrayObject* numpy_array=PyArray_GETCONTIGUOUS((PyArrayObject*)obj_ptr);
    unsigned int array_dim = PyArray_NDIM(numpy_array);
    //require matching dimensions
    if(array_dim!=Dim)
      return(NULL);
    NPY_TYPES type = (NPY_TYPES) PyArray_DESCR(numpy_array)->type_num;
    //require a sane type
    switch(type){
      case NPY_BOOL:
      case NPY_INT8:
      case NPY_INT16:
      case NPY_INT32:
      case NPY_INT64:
      case NPY_UINT8:
      case NPY_UINT16:
      case NPY_UINT32:
      case NPY_UINT64:
      case NPY_FLOAT32:
      case NPY_FLOAT64:
        break;
      default:
        return(NULL);
    }
    
    return(obj_ptr);
  }
  
  static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data){
    PyArrayObject* numpy_array=PyArray_GETCONTIGUOUS((PyArrayObject*)obj_ptr);
    // get numpy array shape and create marray object
#ifdef NPY_1_7_API_VERSION
    npy_intp* array_shape = PyArray_SHAPE(numpy_array);
#else
    npy_intp* array_shape = PyArray_DIMS(numpy_array);
#endif
    std::vector<size_t> dimensions;
    
    unsigned int array_dim = PyArray_NDIM(numpy_array);
    assert(Dim == array_dim && "Non-matching array dimensions.");
    
    for(unsigned int i = 0; i < Dim; i++)
      dimensions.push_back(array_shape[i]);
    
    void* storage=((boost::python::converter::rvalue_from_python_storage<marray<double,Dim>>*)data)->storage.bytes;
    new (storage)marray<double,Dim>;
    data->convertible = storage;
    marray<double,Dim>* oarray=(marray<double,Dim>*)storage;
    
    oarray->resize(dimensions);
    auto it = oarray->begin();
    
    // create numpy iterator
    NpyIter* iter = NpyIter_New(numpy_array, NPY_ITER_READONLY|
                                NPY_ITER_EXTERNAL_LOOP|
                                NPY_ITER_REFS_OK,
                                NPY_KEEPORDER, NPY_NO_CASTING,
                                NULL);
    
    NpyIter_IterNextFunc* iternext = NpyIter_GetIterNext(iter, NULL);
    char** dataptr = NpyIter_GetDataPtrArray(iter);
    npy_intp* strideptr = NpyIter_GetInnerStrideArray(iter);
    npy_intp* sizeptr = NpyIter_GetInnerLoopSizePtr(iter);
    npy_intp iop, nop = NpyIter_GetNOp(iter);
    
    NPY_TYPES type = (NPY_TYPES) PyArray_DESCR(numpy_array)->type_num;
    
    do{
      char* data = *dataptr;
      npy_intp count = *sizeptr;
      npy_intp stride = *strideptr;
      
      while (count--){
        for (iop = 0; iop < nop; ++iop, data+=stride){
          switch(type){
            case NPY_BOOL: *it++ = *(bool*)data; break;
            case NPY_INT8: *it++ = *(int8_t*)data; break;
            case NPY_INT16: *it++ = *(int16_t*)data; break;
            case NPY_INT32: *it++ = *(int32_t*)data; break;
            case NPY_INT64: *it++ = *(int64_t*)data; break;
            case NPY_UINT8: *it++ = *(uint8_t*)data; break;
            case NPY_UINT16: *it++ = *(uint16_t*)data; break;
            case NPY_UINT32: *it++ = *(uint32_t*)data; break;
            case NPY_UINT64: *it++ = *(uint64_t*)data; break;
            case NPY_FLOAT32: *it++ = *(float*)data; break;
            case NPY_FLOAT64: *it++ = *(double*)data; break;
            default:
              throw std::runtime_error("Unsupported array data type");
          }
        }
      }
    } while(iternext(iter));
    
    NpyIter_Deallocate(iter);
  }
};

// nuSQUIDS wrap functions
/*
static void wrap_WriteStateHDF5(nuSQUIDS* nusq, std::string path){
  nusq->WriteStateHDF5(path);
}

*/
static void wrap_ReadStateHDF5(nuSQUIDS* nusq,std::string hdf5_filename,std::string group = "/", std::string cross_section_grp_loc = ""){
  nusq->ReadStateHDF5(hdf5_filename,group,cross_section_grp_loc);
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

static void wrap_nusqatm_Set_GSL_STEP(nuSQUIDSAtm<>* nusq, GSL_STEP_FUNCTIONS step_enum){
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

// overloaded function macro template creator //
#define MAKE_OVERLOAD_TEMPLATE(name, fname, min_args, max_args) \
template<typename T> \
struct name : \
public boost::python::detail::overloads_common<name<T>>{ \
	BOOST_PYTHON_GEN_MEM_FUNCTION(fname, non_void_return_type, \
		max_args, BOOST_PP_SUB_D(1, max_args, min_args), return) \
	typedef non_void_return_type void_return_type; \
	BOOST_PYTHON_OVERLOAD_CONSTRUCTORS(name, max_args + 1, \
		BOOST_PP_SUB_D(1, max_args, min_args)) \
};

// nuSQUIDS-like overloads factories
MAKE_OVERLOAD_TEMPLATE(WriteStateHDF5Overload,WriteStateHDF5,1,5)
MAKE_OVERLOAD_TEMPLATE(ReadStateHDF5Overload,ReadStateHDF5,1,3)
MAKE_OVERLOAD_TEMPLATE(SetInitialStateH5Overload,Set_initial_state,1,2)

// nuSQUIDSpy module definitions
template<typename BaseType, typename = typename std::enable_if<std::is_base_of<nuSQUIDS,BaseType>::value>::type >
  struct RegisterBasicNuSQuIDSPythonBindings {
    const std::string class_label;
    std::shared_ptr<class_<BaseType, boost::noncopyable, std::shared_ptr<BaseType>>> class_object;
    RegisterBasicNuSQuIDSPythonBindings(std::string class_label):class_label(class_label){
      class_object = std::make_shared<class_<BaseType, boost::noncopyable, std::shared_ptr<BaseType>>>(class_label.c_str(), init<>());

      class_object->def(init<marray<double,1>,unsigned int>(args("E_vector","numneu")));
      class_object->def(init<marray<double,1>,unsigned int,NeutrinoType>(args("E_vector","numneu","NT")));
      class_object->def(init<marray<double,1>,unsigned int,NeutrinoType,bool>(args("E_vector","numneu","NT","iinteraction")));
      class_object->def(init<marray<double,1>,unsigned int,NeutrinoType,bool,std::shared_ptr<NeutrinoCrossSections>>(args("E_vector","numneu","NT","iinteraction","ncs")));
      class_object->def(init<std::string>(args("filename")));
      class_object->def(init<std::string, std::string>(args("filename","root group name")));
      class_object->def(init<std::string, std::string, std::shared_ptr<nusquids::nuSQUIDS::InteractionStructure>>(args("filename","root group name","interaction structure")));
      class_object->def(init<unsigned int,NeutrinoType>(args("numneu","NT")));
      class_object->def("Set_initial_state",(void(BaseType::*)(const marray<double,1>&, Basis))&BaseType::Set_initial_state,SetInitialStateH5Overload<BaseType>());
      class_object->def("Set_initial_state",(void(BaseType::*)(const marray<double,2>&, Basis))&BaseType::Set_initial_state,SetInitialStateH5Overload<BaseType>());
      class_object->def("Set_initial_state",(void(BaseType::*)(const marray<double,3>&, Basis))&BaseType::Set_initial_state,SetInitialStateH5Overload<BaseType>());
      class_object->def("Set_Body",&BaseType::Set_Body, bp::arg("Body"));
      class_object->def("Set_Track",&BaseType::Set_Track, bp::arg("Track"));
      class_object->def("Set_E",&BaseType::Set_E, bp::arg("NeutrinoEnergy"));
      class_object->def("EvolveState",&BaseType::EvolveState);
      class_object->def("GetERange",&BaseType::GetERange);
      class_object->def("WriteStateHDF5",&BaseType::WriteStateHDF5,
          WriteStateHDF5Overload<BaseType>(args("hdf5_filename","group"," save_cross_sections","cross_section_grp_loc","overwrite"),
            "Writes the current object into an HDF5 file."));
      class_object->def("ReadStateHDF5",&BaseType::ReadStateHDF5,
          ReadStateHDF5Overload<BaseType>(args("hdf5_filename","group","cross_section_grp_loc"),
            "Reads an HDF5 file and loads the contents into the current object."));
      class_object->def("GetNumNeu",&BaseType::GetNumNeu);
      class_object->def("EvalMass",(double(BaseType::*)(unsigned int) const)&BaseType::EvalMass);
      class_object->def("EvalMass",(double(BaseType::*)(unsigned int,double,unsigned int,double, std::vector<bool>&) const)&BaseType::EvalMass);
      class_object->def("EvalFlavor",(double(BaseType::*)(unsigned int) const)&BaseType::EvalFlavor);
      class_object->def("EvalMass",(double(BaseType::*)(unsigned int,double,unsigned int) const)&BaseType::EvalMass);
      class_object->def("EvalFlavor",(double(BaseType::*)(unsigned int,double,unsigned int) const)&BaseType::EvalFlavor);
      class_object->def("EvalFlavor",(double(BaseType::*)(unsigned int,double,unsigned int,double, std::vector<bool>&) const)&BaseType::EvalFlavor);
      class_object->def("EvalMassAtNode",(double(BaseType::*)(unsigned int,unsigned int,unsigned int) const)&BaseType::EvalMassAtNode);
      class_object->def("EvalFlavorAtNode",(double(BaseType::*)(unsigned int,unsigned int,unsigned int) const)&BaseType::EvalFlavorAtNode);
      class_object->def("GetHamiltonian",&BaseType::GetHamiltonian);
      //class_object->def("GetState",&BaseType::GetState);
      class_object->def("Set_h_min",&BaseType::Set_h_min);
      class_object->def("Set_h_max",&BaseType::Set_h_max);
      class_object->def("Set_h",&BaseType::Set_h);
      class_object->def("Set_rel_error",&BaseType::Set_rel_error);
      class_object->def("Set_abs_error",&BaseType::Set_abs_error);
      class_object->def("Set_AdaptiveStep",&BaseType::Set_AdaptiveStep);
      class_object->def("Set_GSL_step",wrap_Set_GSL_STEP);
      class_object->def("Set_TauRegeneration",&BaseType::Set_TauRegeneration);
      class_object->def("Set_GlashowResonance",&BaseType::Set_GlashowResonance);
      class_object->def("Set_IncludeOscillations",&BaseType::Set_IncludeOscillations);
      class_object->def("Set_AllowConstantDensityOscillationOnlyEvolution",&BaseType::Set_AllowConstantDensityOscillationOnlyEvolution);
      class_object->def("Set_PositivityConstrain",&BaseType::Set_PositivityConstrain);
      class_object->def("Set_PositivityConstrainStep",&BaseType::Set_PositivityConstrainStep);
      class_object->def("Set_ProgressBar",&BaseType::Set_ProgressBar);
      class_object->def("Set_MixingParametersToDefault",&BaseType::Set_MixingParametersToDefault);
      class_object->def("Set_Basis",&BaseType::Set_Basis);
      class_object->def("Set_MixingAngle",&BaseType::Set_MixingAngle);
      class_object->def("Get_MixingAngle",&BaseType::Get_MixingAngle);
      class_object->def("Set_CPPhase",&BaseType::Set_CPPhase);
      class_object->def("Get_CPPhase",&BaseType::Get_CPPhase);
      class_object->def("Set_SquareMassDifference",&BaseType::Set_SquareMassDifference);
      class_object->def("Get_SquareMassDifference",&BaseType::Get_SquareMassDifference);
      class_object->def("GetERange",&BaseType::GetERange);
      class_object->def("GetTrack",&BaseType::GetTrack);
      class_object->def("GetBody",&BaseType::GetBody);
      class_object->def("GetNumE",&BaseType::GetNumE);
      class_object->def("GetNumRho",&BaseType::GetNumRho);
      class_object->def("GetUseInteractions",&BaseType::GetUseInteractions);
      class_object->def("GetUseOscillations",&BaseType::GetUseOscillations);
      class_object->def("InitializeInteractions",&BaseType::InitializeInteractions);
      class_object->def("GetInteractionStructure",(std::shared_ptr<nusquids::nuSQUIDS::InteractionStructure>(BaseType::*)())&BaseType::GetInteractionStructure);
      class_object->def("GetHamiltonian",&BaseType::GetHamiltonian);
      class_object->def("GetTransformationMatrix",&BaseType::GetTransformationMatrix);
      class_object->def("GetNeutrinoCrossSections",&BaseType::GetNeutrinoCrossSections);
      class_object->def("SetNeutrinoCrossSections",&BaseType::SetNeutrinoCrossSections);
      class_object->def("Set_Debug",&BaseType::Set_Debug);
      class_object->def("Set_IncludeOscillations",&BaseType::Set_IncludeOscillations);
      class_object->def("Set_GlashowResonance",&BaseType::Set_GlashowResonance);
    }
    std::shared_ptr<class_<BaseType, boost::noncopyable, std::shared_ptr<BaseType>>> GetClassObject() {
      return class_object;
    }
};

// nuSQUIDSAtm-like overloads factories
MAKE_OVERLOAD_TEMPLATE(nuSQUIDSAtm_EvalFlavor_overload,EvalFlavor,3,5)
MAKE_OVERLOAD_TEMPLATE(nuSQUIDSAtm_Set_initial_state,Set_initial_state,1,2)

// registration for atmospheric template
template<typename BaseType, typename = typename std::enable_if<std::is_base_of<nuSQUIDS,BaseType>::value>::type >
  struct RegisterBasicAtmNuSQuIDSPythonBindings {
    const std::string class_label;
    std::shared_ptr<class_<nuSQUIDSAtm<BaseType>, boost::noncopyable, std::shared_ptr<nuSQUIDSAtm<BaseType>>>> class_object;
    RegisterBasicAtmNuSQuIDSPythonBindings(std::string class_label){
      class_object = std::make_shared<class_<nuSQUIDSAtm<BaseType>, boost::noncopyable, std::shared_ptr<nuSQUIDSAtm<BaseType>>>>(class_label.c_str(), no_init);

      class_object->def(init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType>(args("CosZenith_vector","E_vector","numneu","NT")));
      class_object->def(init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType,bool>(args("CosZenith_vector","E_vector","numneu","NT","iinteraction")));
      class_object->def(init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType,bool,std::shared_ptr<NeutrinoCrossSections>>(args("CosZenith_vector","E_vector","numneu","NT","iinteraction","ncs")));
      class_object->def(init<std::string>(args("filename")));
      class_object->def("EvolveState",&nuSQUIDSAtm<BaseType>::EvolveState);
      class_object->def("Set_TauRegeneration",&nuSQUIDSAtm<BaseType>::Set_TauRegeneration);
      class_object->def("EvalFlavor",(double(nuSQUIDSAtm<BaseType>::*)(unsigned int,double,double,unsigned int,bool) const)&nuSQUIDSAtm<BaseType>::EvalFlavor,
          nuSQUIDSAtm_EvalFlavor_overload<nuSQUIDSAtm<BaseType>>(args("Flavor","cos(theta)","Neutrino Energy","NeuType","BoolToRandomzeProdutionHeight"),
            "nuSQuIDSAtm evaluate flux.."));
      class_object->def("Set_EvalThreads",&nuSQUIDSAtm<BaseType>::Set_EvalThreads);
      class_object->def("Get_EvalThreads",&nuSQUIDSAtm<BaseType>::Get_EvalThreads);
      class_object->def("Set_EarthModel",&nuSQUIDSAtm<BaseType>::Set_EarthModel);
      class_object->def("WriteStateHDF5",&nuSQUIDSAtm<BaseType>::WriteStateHDF5);
      class_object->def("ReadStateHDF5",&nuSQUIDSAtm<BaseType>::ReadStateHDF5);
      class_object->def("Set_MixingAngle",&nuSQUIDSAtm<BaseType>::Set_MixingAngle);
      class_object->def("Get_MixingAngle",&nuSQUIDSAtm<BaseType>::Get_MixingAngle);
      class_object->def("Set_CPPhase",&nuSQUIDSAtm<BaseType>::Set_CPPhase);
      class_object->def("Get_CPPhase",&nuSQUIDSAtm<BaseType>::Get_CPPhase);
      class_object->def("Set_SquareMassDifference",&nuSQUIDSAtm<BaseType>::Set_SquareMassDifference);
      class_object->def("Get_SquareMassDifference",&nuSQUIDSAtm<BaseType>::Get_SquareMassDifference);
      class_object->def("Set_h",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_h);
      class_object->def("Set_h",(void(nuSQUIDSAtm<BaseType>::*)(double,unsigned int))&nuSQUIDSAtm<BaseType>::Set_h);
      class_object->def("Set_h_max",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_h_max);
      class_object->def("Set_h_max",(void(nuSQUIDSAtm<BaseType>::*)(double,unsigned int))&nuSQUIDSAtm<BaseType>::Set_h_max);
      class_object->def("Set_h_min",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_h_min);
      class_object->def("Set_h_min",(void(nuSQUIDSAtm<BaseType>::*)(double,unsigned int))&nuSQUIDSAtm<BaseType>::Set_h_min);
      class_object->def("Set_ProgressBar",&nuSQUIDSAtm<BaseType>::Set_ProgressBar);
      class_object->def("Set_MixingParametersToDefault",&nuSQUIDSAtm<BaseType>::Set_MixingParametersToDefault);
      class_object->def("Set_GSL_step",wrap_nusqatm_Set_GSL_STEP);
      class_object->def("Set_rel_error",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_rel_error);
      class_object->def("Set_rel_error",(void(nuSQUIDSAtm<BaseType>::*)(double, unsigned int))&nuSQUIDSAtm<BaseType>::Set_rel_error);
      class_object->def("Set_abs_error",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_abs_error);
      class_object->def("Set_abs_error",(void(nuSQUIDSAtm<BaseType>::*)(double, unsigned int))&nuSQUIDSAtm<BaseType>::Set_abs_error);
      class_object->def("GetNumE",&nuSQUIDSAtm<BaseType>::GetNumE);
      class_object->def("GetNumCos",&nuSQUIDSAtm<BaseType>::GetNumCos);
      class_object->def("GetNumNeu",&nuSQUIDSAtm<BaseType>::GetNumNeu);
      class_object->def("GetNumRho",&nuSQUIDSAtm<BaseType>::GetNumRho);
      class_object->def("GetnuSQuIDS",(std::vector<BaseType>&(nuSQUIDSAtm<BaseType>::*)())&nuSQUIDSAtm<BaseType>::GetnuSQuIDS,boost::python::return_internal_reference<>());
      class_object->def("GetnuSQuIDS",(BaseType&(nuSQUIDSAtm<BaseType>::*)(unsigned int))&nuSQUIDSAtm<BaseType>::GetnuSQuIDS,boost::python::return_internal_reference<>());
      class_object->def("Set_initial_state",(void(nuSQUIDSAtm<BaseType>::*)(const marray<double,3>&, Basis))&nuSQUIDSAtm<BaseType>::Set_initial_state,nuSQUIDSAtm_Set_initial_state<nuSQUIDSAtm<BaseType>>());
      class_object->def("Set_initial_state",(void(nuSQUIDSAtm<BaseType>::*)(const marray<double,4>&, Basis))&nuSQUIDSAtm<BaseType>::Set_initial_state,nuSQUIDSAtm_Set_initial_state<nuSQUIDSAtm<BaseType>>());
      class_object->def("GetERange",&nuSQUIDSAtm<BaseType>::GetERange);
      class_object->def("GetCosthRange",&nuSQUIDSAtm<BaseType>::GetCosthRange);
      class_object->def("Set_IncludeOscillations",&nuSQUIDSAtm<BaseType>::Set_IncludeOscillations);
      class_object->def("Set_GlashowResonance",&nuSQUIDSAtm<BaseType>::Set_GlashowResonance);
      class_object->def("Set_TauRegeneration",&nuSQUIDSAtm<BaseType>::Set_TauRegeneration);
      class_object->def("Set_AllowConstantDensityOscillationOnlyEvolution",&nuSQUIDSAtm<BaseType>::Set_AllowConstantDensityOscillationOnlyEvolution);
      class_object->def("Set_PositivyConstrain",&nuSQUIDSAtm<BaseType>::Set_PositivityConstrain);
      class_object->def("Set_PositivyConstrainStep",&nuSQUIDSAtm<BaseType>::Set_PositivityConstrainStep);
      class_object->def("Get_EvalThreads",&nuSQUIDSAtm<BaseType>::Get_EvalThreads);
      class_object->def("Set_EvalThreads",&nuSQUIDSAtm<BaseType>::Set_EvalThreads);
      class_object->def("Set_EarthModel",&nuSQUIDSAtm<BaseType>::Set_EarthModel);
      class_object->def("SetNeutrinoCrossSections",&nuSQUIDSAtm<BaseType>::SetNeutrinoCrossSections);
      class_object->def("GetNeutrinoCrossSections",&nuSQUIDSAtm<BaseType>::GetNeutrinoCrossSections);
    }
    std::shared_ptr<class_<nuSQUIDSAtm<BaseType>, boost::noncopyable, std::shared_ptr<nuSQUIDSAtm<BaseType>>>> GetClassObject() {
      return class_object;
    }
};

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

  enum_<Basis>("Basis")
    .value("mass",mass)
    .value("flavor",flavor)
    .value("interaction",interaction)
  ;

  enum_<NeutrinoCrossSections::NeutrinoFlavor>("NeutrinoCrossSections_NeutrinoFlavor")
    .value("electron",NeutrinoCrossSections::NeutrinoFlavor::electron)
    .value("muon",NeutrinoCrossSections::NeutrinoFlavor::muon)
    .value("tau",NeutrinoCrossSections::NeutrinoFlavor::tau)
    .value("sterile",NeutrinoCrossSections::NeutrinoFlavor::sterile)
  ;

  enum_<NeutrinoCrossSections::NeutrinoType>("NeutrinoCrossSections_NeutrinoType")
    .value("neutrino",NeutrinoCrossSections::NeutrinoType::neutrino)
    .value("antineutrino",NeutrinoCrossSections::NeutrinoType::antineutrino)
  ;

  enum_<NeutrinoCrossSections::Current>("NeutrinoCrossSections_Current")
    .value("CC",NeutrinoCrossSections::Current::CC)
    .value("NC",NeutrinoCrossSections::Current::NC)
    .value("GR",NeutrinoCrossSections::Current::GR)
  ;

  class_<squids::SU_vector, boost::noncopyable,std::shared_ptr<squids::SU_vector> >("SU_vector")
    .def(init< std::vector<double> >())
    .def(init<unsigned int>())
    //.def("Rotate",&squids::SU_vector::Rotate)
    .def("Dim",&squids::SU_vector::Dim)
    .def("GetComponents",&squids::SU_vector::GetComponents)
  ;

  enum_<NeutrinoType>("NeutrinoType")
    .value("neutrino",neutrino)
    .value("antineutrino",antineutrino)
    .value("both",both)
  ;

  RegisterBasicNuSQuIDSPythonBindings<nuSQUIDS>("nuSQUIDS");
  RegisterBasicAtmNuSQuIDSPythonBindings<nuSQUIDS>("nuSQUIDSAtm");

  class_<NeutrinoCrossSections, std::shared_ptr<NeutrinoCrossSections>, boost::noncopyable >("NeutrinoCrossSections", no_init);

  class_<NullCrossSections, bases<NeutrinoCrossSections>, std::shared_ptr<NullCrossSections>, boost::noncopyable >("NullCrossSections")
    .def("TotalCrossSection",&NullCrossSections::TotalCrossSection)
    .def("SingleDifferentialCrossSection",&NullCrossSections::SingleDifferentialCrossSection)
  ;

  class_<GlashowResonanceCrossSection, bases<NeutrinoCrossSections>, std::shared_ptr<GlashowResonanceCrossSection>, boost::noncopyable >("GlashowResonanceCrossSection")
    .def("TotalCrossSection",&GlashowResonanceCrossSection::TotalCrossSection)
    .def("SingleDifferentialCrossSection",&GlashowResonanceCrossSection::SingleDifferentialCrossSection)
    .def("WDecayBranchingFraction",&GlashowResonanceCrossSection::WDecayBranchingFraction)
  ;

  class_<NeutrinoDISCrossSectionsFromTables, bases<NeutrinoCrossSections>, std::shared_ptr<NeutrinoDISCrossSectionsFromTables>, boost::noncopyable>("NeutrinoDISCrossSectionsFromTables")
    .def(init<std::string>())
    .def("TotalCrossSection",&NeutrinoDISCrossSectionsFromTables::TotalCrossSection)
    .def("SingleDifferentialCrossSection",&NeutrinoDISCrossSectionsFromTables::SingleDifferentialCrossSection)
    .def("GetNumE",&NeutrinoDISCrossSectionsFromTables::GetNumE)
    .def("GetEmin",&NeutrinoDISCrossSectionsFromTables::GetEmax)
    .def("IsInit",&NeutrinoDISCrossSectionsFromTables::IsInit)
    .def("WriteHDF",&NeutrinoDISCrossSectionsFromTables::WriteHDF)
    .def("WriteText",&NeutrinoDISCrossSectionsFromTables::WriteText)
  ;

  class_<squids::Const, boost::noncopyable>("Const")
    .def_readonly("PeV",&squids::Const::PeV)
    .def_readonly("TeV",&squids::Const::TeV)
    .def_readonly("GeV",&squids::Const::GeV)
    .def_readonly("MeV",&squids::Const::MeV)
    .def_readonly("keV",&squids::Const::keV)
    .def_readonly("eV",&squids::Const::eV)
    .def_readonly("kg",&squids::Const::kg)
    .def_readonly("gr",&squids::Const::gr)
    .def_readonly("meter",&squids::Const::meter)
    .def_readonly("cm",&squids::Const::cm)
    .def_readonly("km",&squids::Const::km)
    .def_readonly("fermi",&squids::Const::fermi)
    .def_readonly("angstrom",&squids::Const::angstrom)
    .def_readonly("AU",&squids::Const::AU)
    .def_readonly("parsec",&squids::Const::parsec)
    .def_readonly("pb",&squids::Const::picobarn)
    .def_readonly("fb",&squids::Const::femtobarn)
    .def_readonly("sec",&squids::Const::sec)
    .def_readonly("hour",&squids::Const::hour)
    .def_readonly("day",&squids::Const::day)
    .def_readonly("year",&squids::Const::year)
  ;

  {
    scope outer
    = class_<Body, std::shared_ptr<Body>, boost::noncopyable >("Body", no_init)
    .def("density",&Body::density)
    .def("ye",&Body::ye)
    ;

    class_<Body::Track, std::shared_ptr<Body::Track>, boost::noncopyable >("Track", no_init)
    .def("GetInitialX",&Body::Track::GetInitialX)
    .def("GetFinalX",&Body::Track::GetFinalX)
    .def("GetX",&Body::Track::GetX)
    .def("SetX",&Body::Track::SetX)
    .def("ReverseTrack",&Body::Track::ReverseTrack)
    ;
  }

  {
    scope outer
    = class_<Vacuum, bases<Body>, std::shared_ptr<Vacuum> >("Vacuum")
    .def("density",&Vacuum::density)
    .def("ye",&Vacuum::ye)
    ;

    class_<Vacuum::Track, std::shared_ptr<Vacuum::Track> >("Track", init<double>())
    .def(init<double,double>())
    .def("GetInitialX",&Vacuum::Track::GetInitialX)
    .def("GetFinalX",&Vacuum::Track::GetFinalX)
    .def("GetX",&Vacuum::Track::GetX)
    .def("SetX",&Vacuum::Track::SetX)
    .def("ReverseTrack",&Vacuum::Track::ReverseTrack)
    ;

    implicitly_convertible< std::shared_ptr<Vacuum>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<Vacuum::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<ConstantDensity, bases<Body>, std::shared_ptr<ConstantDensity> >("ConstantDensity", init<double,double>())
    .def("density",&ConstantDensity::density)
    .def("ye",&ConstantDensity::ye)
    ;

    class_<ConstantDensity::Track, std::shared_ptr<ConstantDensity::Track> >("Track", init<double>())
    .def(init<double,double>())
    .def("GetInitialX",&ConstantDensity::Track::GetInitialX)
    .def("GetFinalX",&ConstantDensity::Track::GetFinalX)
    .def("GetX",&ConstantDensity::Track::GetX)
    .def("SetX",&ConstantDensity::Track::SetX)
    .def("ReverseTrack",&ConstantDensity::Track::ReverseTrack)
    ;

    implicitly_convertible< std::shared_ptr<ConstantDensity>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<ConstantDensity::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<VariableDensity, bases<Body>, std::shared_ptr<VariableDensity> >("VariableDensity", init< std::vector<double>,std::vector<double>,std::vector<double> >())
    ;

    class_<VariableDensity::Track, std::shared_ptr<VariableDensity::Track> >("Track", init<double>())
    .def(init<double,double>())
    .def("GetInitialX",&VariableDensity::Track::GetInitialX)
    .def("GetFinalX",&VariableDensity::Track::GetFinalX)
    .def("GetX",&VariableDensity::Track::GetX)
    .def("SetX",&VariableDensity::Track::SetX)
    .def("ReverseTrack",&VariableDensity::Track::ReverseTrack)
    ;

    implicitly_convertible< std::shared_ptr<VariableDensity>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<VariableDensity::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<Earth, bases<Body>, std::shared_ptr<Earth> >("Earth")
    .def(init<std::string>())
    ;

    class_<Earth::Track, std::shared_ptr<Earth::Track> >("Track", init<double>())
    .def(init<double,double,double>())
    .def("GetInitialX",&Earth::Track::GetInitialX)
    .def("GetFinalX",&Earth::Track::GetFinalX)
    .def("GetX",&Earth::Track::GetX)
    .def("SetX",&Earth::Track::SetX)
    .def("ReverseTrack",&Earth::Track::ReverseTrack)
    ;

    implicitly_convertible< std::shared_ptr<Earth>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<Earth::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<Sun, bases<Body>, std::shared_ptr<Sun> >("Sun")
    ;

    class_<Sun::Track, std::shared_ptr<Sun::Track> >("Track", init<double>())
    .def(init<double,double>())
    .def("GetInitialX",&Sun::Track::GetInitialX)
    .def("GetFinalX",&Sun::Track::GetFinalX)
    .def("GetX",&Sun::Track::GetX)
    .def("SetX",&Sun::Track::SetX)
    .def("ReverseTrack",&Sun::Track::ReverseTrack)
    ;

    implicitly_convertible< std::shared_ptr<Sun>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<Sun::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<SunASnu, bases<Body>, std::shared_ptr<SunASnu> >("SunASnu")
    ;

    class_<SunASnu::Track, std::shared_ptr<SunASnu::Track> >("Track", init<double>())
    .def(init<double,double>())
    .def("GetInitialX",&SunASnu::Track::GetInitialX)
    .def("GetFinalX",&SunASnu::Track::GetFinalX)
    .def("GetX",&SunASnu::Track::GetX)
    .def("SetX",&SunASnu::Track::SetX)
    .def("ReverseTrack",&SunASnu::Track::ReverseTrack)
    ;

    implicitly_convertible< std::shared_ptr<SunASnu>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<SunASnu::Track>, std::shared_ptr<Body::Track> >();
  }

  {
    scope outer
    = class_<EarthAtm, bases<Body>, boost::noncopyable, std::shared_ptr<EarthAtm> >("EarthAtm")
    .def(init<std::string>())
    ;

    class_<EarthAtm::Track, std::shared_ptr<EarthAtm::Track> >("Track", init<double>())
    .def("GetInitialX",&EarthAtm::Track::GetInitialX)
    .def("GetFinalX",&EarthAtm::Track::GetFinalX)
    .def("GetX",&EarthAtm::Track::GetX)
    .def("SetX",&EarthAtm::Track::SetX)
    .def("ReverseTrack",&EarthAtm::Track::ReverseTrack)
    ;

    implicitly_convertible< std::shared_ptr<EarthAtm>, std::shared_ptr<Body> >();
    implicitly_convertible< std::shared_ptr<EarthAtm::Track>, std::shared_ptr<Body::Track> >();
  }

  // python cotainer to vector<double> convertion
  using namespace scitbx::boost_python::container_conversions;
  from_python_sequence< std::vector<double>, variable_capacity_policy >();
  //from_python_sequence< std::vector< std::vector<double> >, variable_capacity_policy >();
  to_python_converter< std::vector<double, class std::allocator<double> >, VecToList<double> > ();
  // register marray converters
  to_python_converter< marray<double,1> , marray_to_numpyarray<1> >();
  to_python_converter< marray<double,2> , marray_to_numpyarray<2> >();
  to_python_converter< marray<double,3> , marray_to_numpyarray<3> >();
  to_python_converter< marray<double,4> , marray_to_numpyarray<4> >();

  marray_from_python<1>();
  marray_from_python<2>();
  marray_from_python<3>();
  marray_from_python<4>();
}
