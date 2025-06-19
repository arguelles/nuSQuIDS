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

#ifndef NUSQUIDS_PY_H
#define NUSQUIDS_PY_H

#if __cplusplus < 201103L
#error C++11 compiler required. Update your compiler and use the flag -std=c++11
#endif

#define H5Gopen_vers 2
#define H5Gcreate_vers 2
#define H5Eset_auto_vers 2
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/resources.h>
#include <nuSQuIDS/marray.h>

#include <numpy/ndarrayobject.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

using namespace nusquids;

template<class T>
std::string PrintObject(const T& rObject)
{
    std::stringstream ss;
    ss << rObject;
    return ss.str();
}

namespace pybind11 { namespace detail {
	template <unsigned int Dim>
	struct type_caster<marray<double, Dim>> {
	private:
	    using T = marray<double, Dim>;
	public:
            PYBIND11_TYPE_CASTER(T, _("marray<double,Dim>"));

	    // Python -> C++
	    bool load(handle src, bool) {
		if (!py::isinstance<py::array>(src))
		    return false;

		auto array = py::array::ensure(src);
		if (!array || array.ndim() != Dim || !py::isinstance<py::array_t<double>>(array))
		    return false;

		std::array<size_t, Dim> shape;
		for (size_t i = 0; i < Dim; ++i)
		    shape[i] = static_cast<size_t>(array.shape(i));

		double* ptr = static_cast<double*>(array.mutable_data());
		if (!ptr)
		    return false;

		value = T();
		value.resize(shape);
		std::memcpy(value.get_data(), array.data(), sizeof(double) * array.size());
		return true;
	    }

	    // C++ -> Python
	    static handle cast(const marray<double, Dim>& arr, return_value_policy, handle parent) {
		std::vector<ssize_t> shape(Dim);
		std::vector<ssize_t> strides(Dim);

		ssize_t stride = sizeof(double);
		for (ssize_t i = Dim - 1; i >= 0; --i) {
		    shape[i] = arr.extent(i);
		    strides[i] = stride;
		    stride *= shape[i];
		}

		return py::array(py::buffer_info(
		    const_cast<double*>(arr.get_data()), // assume mutable data
		    sizeof(double),
		    py::format_descriptor<double>::format(),
		    Dim,
		    shape,
		    strides
		)).release();
	    }
	};
}} // namespace pybind11::detail

enum GSL_STEP_FUNCTIONS {
  GSL_STEP_RK2,
  GSL_STEP_RK4,
  GSL_STEP_RKF45,
  GSL_STEP_RKCK,
  GSL_STEP_RK8PD,
  GSL_STEP_MSADAMS
};

template<typename BaseType, typename = typename std::enable_if<std::is_base_of<nuSQUIDS,BaseType>::value>::type >
static void wrap_Set_GSL_STEP(BaseType* nusq, GSL_STEP_FUNCTIONS step_enum){
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
    case GSL_STEP_MSADAMS:
      nusq->Set_GSL_step(gsl_odeiv2_step_msadams);
      break;
  }
}

template<typename BaseType, typename = typename std::enable_if<std::is_base_of<nuSQUIDS,BaseType>::value>::type >
static void wrap_nusqatm_Set_GSL_STEP(nuSQUIDSAtm<BaseType>* nusq, GSL_STEP_FUNCTIONS step_enum){
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
    case GSL_STEP_MSADAMS:
      nusq->Set_GSL_step(gsl_odeiv2_step_msadams);
      break;
  }
}

// nuSQUIDSpy module definitions
template<typename BaseType, typename = typename std::enable_if<std::is_base_of<nuSQUIDS,BaseType>::value>::type >
  struct RegisterBasicNuSQuIDSPythonBindings {
    const std::string class_label;
    std::shared_ptr<py::class_<BaseType, std::shared_ptr<BaseType>>> class_object;
    RegisterBasicNuSQuIDSPythonBindings(py::module_ m, std::string class_label):class_label(class_label){
      class_object = std::make_shared<py::class_<BaseType, std::shared_ptr<BaseType>>>(m,class_label.c_str());

      class_object->def(py::init<>());
      class_object->def(py::init<marray<double,1>,unsigned int>(),py::arg("E_vector"),py::arg("numneu"));
      class_object->def(py::init<marray<double,1>,unsigned int,NeutrinoType>(),py::arg("E_vector"),py::arg("numneu"),py::arg("NT"));
      class_object->def(py::init<marray<double,1>,unsigned int,NeutrinoType,bool>(),py::arg("E_vector"),py::arg("numneu"),py::arg("NT"),py::arg("iinteraction"));
      class_object->def(py::init<marray<double,1>,unsigned int,NeutrinoType,bool,std::shared_ptr<CrossSectionLibrary>>(),py::arg("E_vector"),py::arg("numneu"),py::arg("NT"),py::arg("iinteraction"),py::arg("ncs"));
      class_object->def(py::init<std::string>(),py::arg("filename"));
      class_object->def(py::init<std::string, std::string>(),py::arg("filename"),py::arg("root group name"));
      class_object->def(py::init<std::string, std::string, std::shared_ptr<nusquids::nuSQUIDS::InteractionStructure>>(),py::arg("filename"),py::arg("root group name"),py::arg("interaction structure"));
      class_object->def(py::init<unsigned int,NeutrinoType>(),py::arg("numneu"),py::arg("NT"));
      class_object->def("Set_initial_state",(void(BaseType::*)(const marray<double,1>&, Basis))&BaseType::Set_initial_state, py::arg("ini_state"), py::arg("basis") = Basis::flavor);
      class_object->def("Set_initial_state",(void(BaseType::*)(const marray<double,2>&, Basis))&BaseType::Set_initial_state, py::arg("ini_state"), py::arg("basis") = Basis::flavor);
      class_object->def("Set_initial_state",(void(BaseType::*)(const marray<double,3>&, Basis))&BaseType::Set_initial_state, py::arg("ini_state"), py::arg("basis") = Basis::flavor);
      class_object->def("Set_Body",&BaseType::Set_Body, py::arg("Body"));
      class_object->def("Set_Track",&BaseType::Set_Track, py::arg("Track"));
      class_object->def("Set_E",&BaseType::Set_E, py::arg("NeutrinoEnergy"));
      class_object->def("EvolveState",&BaseType::EvolveState);
      class_object->def("GetERange",&BaseType::GetERange);
      class_object->def("WriteStateHDF5",&BaseType::WriteStateHDF5, py::arg("hdf5_filename"),py::arg("group") = "/",py::arg("save_cross_sections") = true,py::arg("cross_section_grp_loc") = "",py::arg("overwrite") = true,
            "Writes the current object into an HDF5 file.");
      class_object->def("ReadStateHDF5",&BaseType::ReadStateHDF5,py::arg("hdf5_filename"),py::arg("group") = "/",py::arg("cross_section_grp_loc") = "",
            "Reads an HDF5 file and loads the contents into the current object.");
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
      class_object->def("GetState",(const squids::SU_vector&(BaseType::*)(unsigned int))&BaseType::GetState, py::return_value_policy::copy);
      class_object->def("GetState",(const squids::SU_vector&(BaseType::*)(unsigned int, unsigned int))&BaseType::GetState, py::return_value_policy::copy);
      class_object->def("Set_EvolLowPassCutoff", &BaseType::Set_EvolLowPassCutoff);
      class_object->def("Set_EvolLowPassScale", &BaseType::Set_EvolLowPassScale);
      class_object->def("Set_h_min",&BaseType::Set_h_min);
      class_object->def("Set_h_max",&BaseType::Set_h_max);
      class_object->def("Set_h",&BaseType::Set_h);
      class_object->def("Set_rel_error",&BaseType::Set_rel_error);
      class_object->def("Set_abs_error",&BaseType::Set_abs_error);
      class_object->def("Set_AdaptiveStep",&BaseType::Set_AdaptiveStep);
      class_object->def("Set_GSL_step",wrap_Set_GSL_STEP<BaseType>);
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
      class_object->def("Set_NeutrinoSources",&BaseType::Set_NeutrinoSources);
      class_object->def("Get_NeutrinoSources",&BaseType::Get_NeutrinoSources);
    }
    std::shared_ptr<py::class_<BaseType, std::shared_ptr<BaseType>>> GetClassObject() {
      return class_object;
    }
};

// registration for atmospheric template
template<typename BaseType, typename = typename std::enable_if<std::is_base_of<nuSQUIDS,BaseType>::value>::type >
  struct RegisterBasicAtmNuSQuIDSPythonBindings {
    const std::string class_label;
    std::shared_ptr<py::class_<nuSQUIDSAtm<BaseType>, std::shared_ptr<nuSQUIDSAtm<BaseType>>>> class_object;
    RegisterBasicAtmNuSQuIDSPythonBindings(py::module_ m,std::string class_label){
      class_object = std::make_shared<py::class_<nuSQUIDSAtm<BaseType>, std::shared_ptr<nuSQUIDSAtm<BaseType>>>>(m,class_label.c_str());

      class_object->def(py::init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType>(),py::arg("CosZenith_vector"),py::arg("E_vector"),py::arg("numneu"),py::arg("NT"));
      class_object->def(py::init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType,bool>(),py::arg("CosZenith_vector"),py::arg("E_vector"),py::arg("numneu"),py::arg("NT"),py::arg("iinteraction"));
      class_object->def(py::init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType,bool,std::shared_ptr<CrossSectionLibrary>>(),py::arg("CosZenith_vector"),py::arg("E_vector"),py::arg("numneu"),py::arg("NT"),py::arg("iinteraction"),py::arg("ncs"));
      class_object->def(py::init<std::string>(),py::arg("filename"));
      class_object->def("EvolveState",&nuSQUIDSAtm<BaseType>::EvolveState);
      class_object->def("Set_TauRegeneration",&nuSQUIDSAtm<BaseType>::Set_TauRegeneration);
      class_object->def("EvalFlavor",(double(nuSQUIDSAtm<BaseType>::*)(unsigned int,double,double,unsigned int,bool) const)&nuSQUIDSAtm<BaseType>::EvalFlavor, py::arg("Flavor"),py::arg("cos(theta)"),py::arg("Neutrino Energy"),py::arg("NeuType") = 0,py::arg("BoolToRandomzeProdutionHeight") = false,
            "nuSQuIDSAtm evaluate flux.");
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
      class_object->def("Set_GSL_step",wrap_nusqatm_Set_GSL_STEP<BaseType>);
      class_object->def("Set_rel_error",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_rel_error);
      class_object->def("Set_rel_error",(void(nuSQUIDSAtm<BaseType>::*)(double, unsigned int))&nuSQUIDSAtm<BaseType>::Set_rel_error);
      class_object->def("Set_abs_error",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_abs_error);
      class_object->def("Set_abs_error",(void(nuSQUIDSAtm<BaseType>::*)(double, unsigned int))&nuSQUIDSAtm<BaseType>::Set_abs_error);
      class_object->def("Set_EvolLowPassCutoff",&nuSQUIDSAtm<BaseType>::Set_EvolLowPassCutoff);
      class_object->def("Set_EvolLowPassScale",&nuSQUIDSAtm<BaseType>::Set_EvolLowPassScale);
      class_object->def("GetNumE",&nuSQUIDSAtm<BaseType>::GetNumE);
      class_object->def("GetNumCos",&nuSQUIDSAtm<BaseType>::GetNumCos);
      class_object->def("GetNumNeu",&nuSQUIDSAtm<BaseType>::GetNumNeu);
      class_object->def("GetNumRho",&nuSQUIDSAtm<BaseType>::GetNumRho);
      class_object->def("GetnuSQuIDS",(std::vector<BaseType>&(nuSQUIDSAtm<BaseType>::*)())&nuSQUIDSAtm<BaseType>::GetnuSQuIDS,py::return_value_policy::reference_internal);
      class_object->def("GetnuSQuIDS",(BaseType&(nuSQUIDSAtm<BaseType>::*)(unsigned int))&nuSQUIDSAtm<BaseType>::GetnuSQuIDS,py::return_value_policy::reference_internal);
      class_object->def("Set_initial_state",(void(nuSQUIDSAtm<BaseType>::*)(const marray<double,3>&, Basis))&nuSQUIDSAtm<BaseType>::Set_initial_state,py::arg("ini_flux"),py::arg("basis") = Basis::flavor);
      class_object->def("Set_initial_state",(void(nuSQUIDSAtm<BaseType>::*)(const marray<double,4>&, Basis))&nuSQUIDSAtm<BaseType>::Set_initial_state,py::arg("ini_flux"),py::arg("basis") = Basis::flavor);
      class_object->def("GetStates", (marray<double,2>(nuSQUIDSAtm<BaseType>::*)(unsigned int))&nuSQUIDSAtm<BaseType>::GetStates,py::arg("rho") = 0, "Get evolved states of all nodes");
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
    std::shared_ptr<py::class_<nuSQUIDSAtm<BaseType>, std::shared_ptr<nuSQUIDSAtm<BaseType>>>> GetClassObject() {
      return class_object;
    }
};

#endif
