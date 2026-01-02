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
      class_object = std::make_shared<py::class_<BaseType, std::shared_ptr<BaseType>>>(m,class_label.c_str(),
R"doc(nuSQUIDS main class for neutrino evolution calculations.

This class solves the neutrino evolution equations in various environments
(vacuum, constant density, Earth, Sun, etc.) and optionally includes
non-coherent interactions (CC, NC, tau regeneration, Glashow resonance).

There are two main modes of operation:
  - Single energy mode: for calculating oscillation probabilities at a single energy
  - Multiple energy mode: for propagating neutrino fluxes across an energy spectrum

Examples
--------
Single energy mode (oscillation probability):
    >>> import nuSQuIDS as nsq
    >>> units = nsq.Const()
    >>> nus = nsq.nuSQUIDS(3, nsq.NeutrinoType.neutrino)
    >>> nus.Set_Body(nsq.Vacuum())
    >>> nus.Set_Track(nsq.Vacuum.Track(100*units.km))
    >>> nus.Set_E(1.0*units.GeV)
    >>> nus.Set_initial_state([0, 1, 0], nsq.Basis.flavor)  # pure nu_mu
    >>> nus.EvolveState()
    >>> print(f"P(nu_mu->nu_e) = {nus.EvalFlavor(0)}")

Multiple energy mode (flux propagation):
    >>> E_nodes = nsq.logspace(1e9, 1e12, 100)  # 1 GeV to 1 TeV
    >>> nus = nsq.nuSQUIDS(E_nodes, 3, nsq.NeutrinoType.both, True)  # with interactions
    >>> # ... set body, track, initial flux, then evolve
)doc");

      // Constructors with docstrings
      class_object->def(py::init<>(),
        "Default constructor. Creates an uninitialized nuSQUIDS object.");

      class_object->def(py::init<marray<double,1>,unsigned int>(),
        py::arg("E_vector"),py::arg("numneu"),
R"doc(Multiple energy mode constructor (neutrino+antineutrino, no interactions).

Parameters
----------
E_vector : array_like
    Energy nodes in eV. Use logspace() or linspace() to generate.
numneu : int
    Number of neutrino flavors (e.g., 3 for standard oscillations, 4+ for sterile).
)doc");

      class_object->def(py::init<marray<double,1>,unsigned int,NeutrinoType>(),
        py::arg("E_vector"),py::arg("numneu"),py::arg("NT"),
R"doc(Multiple energy mode constructor with neutrino type selection.

Parameters
----------
E_vector : array_like
    Energy nodes in eV.
numneu : int
    Number of neutrino flavors.
NT : NeutrinoType
    neutrino, antineutrino, or both (for simultaneous evolution).
)doc");

      class_object->def(py::init<marray<double,1>,unsigned int,NeutrinoType,bool>(),
        py::arg("E_vector"),py::arg("numneu"),py::arg("NT"),py::arg("iinteraction"),
R"doc(Multiple energy mode constructor with interactions.

Parameters
----------
E_vector : array_like
    Energy nodes in eV.
numneu : int
    Number of neutrino flavors.
NT : NeutrinoType
    neutrino, antineutrino, or both.
iinteraction : bool
    If True, include non-coherent interactions (CC, NC).
    Requires NT=both for tau regeneration.
)doc");

      class_object->def(py::init<marray<double,1>,unsigned int,NeutrinoType,bool,std::shared_ptr<CrossSectionLibrary>>(),
        py::arg("E_vector"),py::arg("numneu"),py::arg("NT"),py::arg("iinteraction"),py::arg("ncs"),
R"doc(Multiple energy mode constructor with custom cross sections.

Parameters
----------
E_vector : array_like
    Energy nodes in eV.
numneu : int
    Number of neutrino flavors.
NT : NeutrinoType
    neutrino, antineutrino, or both.
iinteraction : bool
    If True, include non-coherent interactions.
ncs : CrossSectionLibrary
    Custom neutrino cross section library.
)doc");

      class_object->def(py::init<std::string>(),py::arg("filename"),
R"doc(Construct from HDF5 file.

Parameters
----------
filename : str
    Path to HDF5 file containing a saved nuSQUIDS state.
)doc");

      class_object->def(py::init<std::string, std::string>(),
        py::arg("filename"),py::arg("group"),
R"doc(Construct from HDF5 file with custom group.

Parameters
----------
filename : str
    Path to HDF5 file.
group : str
    HDF5 group path where nuSQUIDS state is stored.
)doc");

      class_object->def(py::init<std::string, std::string, std::shared_ptr<nusquids::nuSQUIDS::InteractionStructure>>(),
        py::arg("filename"),py::arg("group"),py::arg("int_struct"),
R"doc(Construct from HDF5 file with pre-computed interaction structure.

Parameters
----------
filename : str
    Path to HDF5 file.
group : str
    HDF5 group path.
int_struct : InteractionStructure
    Pre-computed interaction structure (avoids recomputing cross sections).
)doc");

      class_object->def(py::init<unsigned int,NeutrinoType>(),
        py::arg("numneu"),py::arg("NT"),
R"doc(Single energy mode constructor.

Parameters
----------
numneu : int
    Number of neutrino flavors.
NT : NeutrinoType
    neutrino or antineutrino (both not supported in single energy mode).

Note
----
Interactions are not available in single energy mode.
Use Set_E() to set the neutrino energy before evolving.
)doc");

      // State setting methods
      class_object->def("Set_initial_state",(void(BaseType::*)(const marray<double,1>&, Basis))&BaseType::Set_initial_state,
        py::arg("ini_state"), py::arg("basis") = Basis::flavor,
R"doc(Set initial state in single energy mode.

Parameters
----------
ini_state : array_like
    1D array of length numneu. In flavor basis: [nu_e, nu_mu, nu_tau, ...].
    In mass basis: [nu_1, nu_2, nu_3, ...].
basis : Basis, optional
    Basis of ini_state: flavor (default) or mass.
)doc");

      class_object->def("Set_initial_state",(void(BaseType::*)(const marray<double,2>&, Basis))&BaseType::Set_initial_state,
        py::arg("ini_state"), py::arg("basis") = Basis::flavor,
R"doc(Set initial state in multiple energy mode (neutrino or antineutrino only).

Parameters
----------
ini_state : array_like
    2D array of shape (n_energies, numneu).
    Row i contains flavor/mass composition at energy node i.
basis : Basis, optional
    Basis of ini_state: flavor (default) or mass.
)doc");

      class_object->def("Set_initial_state",(void(BaseType::*)(const marray<double,3>&, Basis))&BaseType::Set_initial_state,
        py::arg("ini_state"), py::arg("basis") = Basis::flavor,
R"doc(Set initial state in multiple energy mode with neutrino and antineutrino.

Parameters
----------
ini_state : array_like
    3D array of shape (n_energies, 2, numneu).
    Axis 1: 0=neutrino, 1=antineutrino.
basis : Basis, optional
    Basis of ini_state: flavor (default) or mass.
)doc");

      class_object->def("Set_Body",&BaseType::Set_Body, py::arg("body"),
R"doc(Set the body (medium) for neutrino propagation.

Parameters
----------
body : Body
    The propagation medium. Options include:
    - Vacuum(): vacuum propagation
    - ConstantDensity(density, ye): constant matter density
    - Earth(): Earth with PREM density profile
    - EarthAtm(): Earth for atmospheric neutrinos
    - Sun(): Solar density profile

Example
-------
>>> nus.Set_Body(nsq.Earth())
)doc");

      class_object->def("Set_Track",&BaseType::Set_Track, py::arg("track"),
R"doc(Set the trajectory through the body.

Parameters
----------
track : Track
    Trajectory object matching the body type.
    Each body has an associated Track class.

Example
-------
>>> nus.Set_Body(nsq.Vacuum())
>>> nus.Set_Track(nsq.Vacuum.Track(1000*units.km))
)doc");

      class_object->def("Set_E",&BaseType::Set_E, py::arg("E"),
R"doc(Set neutrino energy in single energy mode.

Parameters
----------
E : float
    Neutrino energy in eV. Use units for convenience:
    E = 1.0*units.GeV
)doc");

      class_object->def("EvolveState",&BaseType::EvolveState,
R"doc(Evolve the neutrino state through the body along the track.

This is the main evolution method. Call after setting:
- Body (Set_Body)
- Track (Set_Track)
- Initial state (Set_initial_state)
- Energy (Set_E, single energy mode only)

The evolution uses the GSL ODE solver with adaptive step size.
)doc");

      class_object->def("GetERange",&BaseType::GetERange,
R"doc(Get the energy nodes.

Returns
-------
array
    1D array of energy values in eV.
)doc");

      class_object->def("WriteStateHDF5",&BaseType::WriteStateHDF5,
        py::arg("hdf5_filename"),py::arg("group") = "/",py::arg("save_cross_sections") = true,
        py::arg("cross_section_grp_loc") = "",py::arg("overwrite") = true,
R"doc(Write current state to HDF5 file.

Parameters
----------
hdf5_filename : str
    Output filename.
group : str, optional
    HDF5 group path (default: root "/").
save_cross_sections : bool, optional
    If True, save cross section tables (default: True).
cross_section_grp_loc : str, optional
    Custom location for cross sections.
overwrite : bool, optional
    If True, overwrite existing file (default: True).
)doc");

      class_object->def("ReadStateHDF5",&BaseType::ReadStateHDF5,
        py::arg("hdf5_filename"),py::arg("group") = "/",py::arg("cross_section_grp_loc") = "",
R"doc(Read state from HDF5 file.

Parameters
----------
hdf5_filename : str
    Input filename.
group : str, optional
    HDF5 group path (default: root "/").
cross_section_grp_loc : str, optional
    Custom location for cross sections.
)doc");

      class_object->def("GetNumNeu",&BaseType::GetNumNeu,
R"doc(Get the number of neutrino flavors.

Returns
-------
int
    Number of neutrino flavors.
)doc");

      // Evaluation methods
      class_object->def("EvalMass",(double(BaseType::*)(unsigned int) const)&BaseType::EvalMass,
        py::arg("mass_state"),
R"doc(Evaluate mass eigenstate composition (single energy mode).

Parameters
----------
mass_state : int
    Mass eigenstate index (0=nu_1, 1=nu_2, 2=nu_3, ...).

Returns
-------
float
    The density matrix projection onto the mass eigenstate.
)doc");

      class_object->def("EvalMass",(double(BaseType::*)(unsigned int,double,unsigned int,double, std::vector<bool>&) const)&BaseType::EvalMass,
        py::arg("mass_state"),py::arg("E"),py::arg("rho"),py::arg("scale"),py::arg("avr"));

      class_object->def("EvalFlavor",(double(BaseType::*)(unsigned int) const)&BaseType::EvalFlavor,
        py::arg("flavor"),
R"doc(Evaluate flavor composition (single energy mode).

Parameters
----------
flavor : int
    Flavor index (0=nu_e, 1=nu_mu, 2=nu_tau, 3+=sterile).

Returns
-------
float
    The density matrix projection onto the flavor state.
    For oscillation probability, this equals P(initial -> flavor).
)doc");

      class_object->def("EvalMass",(double(BaseType::*)(unsigned int,double,unsigned int) const)&BaseType::EvalMass,
        py::arg("mass_state"),py::arg("E"),py::arg("rho"),
R"doc(Evaluate mass composition at given energy (multiple energy mode).

Parameters
----------
mass_state : int
    Mass eigenstate index.
E : float
    Energy in eV.
rho : int
    Equation index: 0=neutrino, 1=antineutrino (when NT=both).

Returns
-------
float
    Interpolated mass eigenstate content at energy E.
)doc");

      class_object->def("EvalFlavor",(double(BaseType::*)(unsigned int,double,unsigned int) const)&BaseType::EvalFlavor,
        py::arg("flavor"),py::arg("E"),py::arg("rho") = 0,
R"doc(Evaluate flavor composition at given energy (multiple energy mode).

Parameters
----------
flavor : int
    Flavor index (0=nu_e, 1=nu_mu, 2=nu_tau).
E : float
    Energy in eV.
rho : int, optional
    Equation index: 0=neutrino (default), 1=antineutrino.

Returns
-------
float
    Interpolated flux of the given flavor at energy E.

Example
-------
>>> phi_mu = nus.EvalFlavor(1, 1e10, 0)  # nu_mu flux at 10 GeV
)doc");

      class_object->def("EvalFlavor",(double(BaseType::*)(unsigned int,double,unsigned int,double, std::vector<bool>&) const)&BaseType::EvalFlavor,
        py::arg("flavor"),py::arg("E"),py::arg("rho"),py::arg("scale"),py::arg("avr"),
R"doc(Evaluate flavor with oscillation averaging.

Parameters
----------
flavor : int
    Flavor index.
E : float
    Energy in eV.
rho : int
    Equation index.
scale : float
    Scale for oscillation averaging.
avr : list of bool
    Output: which oscillation scales were averaged.

Returns
-------
float
    Flavor content with fast oscillations averaged out.
)doc");

      class_object->def("EvalMassAtNode",(double(BaseType::*)(unsigned int,unsigned int,unsigned int) const)&BaseType::EvalMassAtNode,
        py::arg("mass_state"),py::arg("ie"),py::arg("rho") = 0,
R"doc(Evaluate mass composition at a specific energy node.

Parameters
----------
mass_state : int
    Mass eigenstate index.
ie : int
    Energy node index.
rho : int, optional
    Equation index (default: 0).

Returns
-------
float
    Mass eigenstate content at the node (no interpolation).
)doc");

      class_object->def("EvalFlavorAtNode",(double(BaseType::*)(unsigned int,unsigned int,unsigned int) const)&BaseType::EvalFlavorAtNode,
        py::arg("flavor"),py::arg("ie"),py::arg("rho") = 0,
R"doc(Evaluate flavor composition at a specific energy node.

Parameters
----------
flavor : int
    Flavor index.
ie : int
    Energy node index.
rho : int, optional
    Equation index (default: 0).

Returns
-------
float
    Flavor content at the node (no interpolation).
)doc");

      class_object->def("GetHamiltonian",&BaseType::GetHamiltonian,
        py::arg("ie"),py::arg("rho") = 0,
R"doc(Get the Hamiltonian at an energy node.

Parameters
----------
ie : int
    Energy node index.
rho : int, optional
    Equation index (default: 0).

Returns
-------
SU_vector
    The full Hamiltonian (H0 + matter potential).
)doc");

      class_object->def("GetState",(const squids::SU_vector&(BaseType::*)(unsigned int))&BaseType::GetState,
        py::return_value_policy::copy, py::arg("ie"),
R"doc(Get the density matrix state at an energy node.

Parameters
----------
ie : int
    Energy node index.

Returns
-------
SU_vector
    The density matrix in SU(N) representation.
)doc");

      class_object->def("GetState",(const squids::SU_vector&(BaseType::*)(unsigned int, unsigned int))&BaseType::GetState,
        py::return_value_policy::copy, py::arg("ie"),py::arg("rho"),
R"doc(Get the density matrix state at an energy node.

Parameters
----------
ie : int
    Energy node index.
rho : int
    Equation index.

Returns
-------
SU_vector
    The density matrix in SU(N) representation.
)doc");

      // Evolution control methods
      class_object->def("Set_EvolLowPassCutoff", &BaseType::Set_EvolLowPassCutoff, py::arg("cutoff"),
R"doc(Set cutoff for state evolution low-pass filter.

Parameters
----------
cutoff : float
    Frequency cutoff value.
)doc");

      class_object->def("Set_EvolLowPassScale", &BaseType::Set_EvolLowPassScale, py::arg("scale"),
R"doc(Set scale for state evolution low-pass filter ramp.

Parameters
----------
scale : float
    Frequency range for linear ramp.
)doc");

      class_object->def("Set_h_min",&BaseType::Set_h_min, py::arg("h_min"),
R"doc(Set minimum step size for ODE solver.

Parameters
----------
h_min : float
    Minimum step size in natural units.
)doc");

      class_object->def("Set_h_max",&BaseType::Set_h_max, py::arg("h_max"),
R"doc(Set maximum step size for ODE solver.

Parameters
----------
h_max : float
    Maximum step size in natural units.
    Use units, e.g., 500*units.km.
)doc");

      class_object->def("Set_h",&BaseType::Set_h, py::arg("h"),
R"doc(Set current step size for ODE solver.

Parameters
----------
h : float
    Step size in natural units.
)doc");

      class_object->def("Set_rel_error",&BaseType::Set_rel_error, py::arg("rel_error"),
R"doc(Set relative error tolerance for ODE solver.

Parameters
----------
rel_error : float
    Relative error tolerance (e.g., 1e-9).
)doc");

      class_object->def("Set_abs_error",&BaseType::Set_abs_error, py::arg("abs_error"),
R"doc(Set absolute error tolerance for ODE solver.

Parameters
----------
abs_error : float
    Absolute error tolerance (e.g., 1e-9).
)doc");

      class_object->def("Set_AdaptiveStep",&BaseType::Set_AdaptiveStep, py::arg("adaptive"),
R"doc(Enable or disable adaptive step size.

Parameters
----------
adaptive : bool
    If True, use adaptive step size (recommended).
)doc");

      class_object->def("Set_GSL_step",wrap_Set_GSL_STEP<BaseType>, py::arg("step_type"),
R"doc(Set the GSL ODE stepping algorithm.

Parameters
----------
step_type : GSL_STEP_FUNCTIONS
    Stepping algorithm:
    - GSL_STEP_RK2: 2nd order Runge-Kutta
    - GSL_STEP_RK4: 4th order Runge-Kutta
    - GSL_STEP_RKF45: Runge-Kutta-Fehlberg (4,5) (default)
    - GSL_STEP_RKCK: Runge-Kutta Cash-Karp (4,5)
    - GSL_STEP_RK8PD: Runge-Kutta Prince-Dormand (8,9)
    - GSL_STEP_MSADAMS: Adams multistep
)doc");

      // Physics toggles
      class_object->def("Set_TauRegeneration",&BaseType::Set_TauRegeneration, py::arg("opt"),
R"doc(Toggle tau regeneration.

Parameters
----------
opt : bool
    If True, include tau regeneration from tau decay.
    Requires interactions enabled and NT=both.
)doc");

      class_object->def("Set_GlashowResonance",&BaseType::Set_GlashowResonance, py::arg("opt"),
R"doc(Toggle Glashow resonance.

Parameters
----------
opt : bool
    If True, include W- production via nu_e_bar + e -> W-.
    Important for PeV-scale electron antineutrinos.
)doc");

      class_object->def("Set_IncludeOscillations",&BaseType::Set_IncludeOscillations, py::arg("opt"),
R"doc(Toggle neutrino oscillations.

Parameters
----------
opt : bool
    If True (default), include oscillations.
    Set False for pure absorption/interaction studies.
)doc");

      class_object->def("Set_AllowConstantDensityOscillationOnlyEvolution",&BaseType::Set_AllowConstantDensityOscillationOnlyEvolution, py::arg("opt"),
R"doc(Toggle fast constant density evolution.

Parameters
----------
opt : bool
    If True, use analytic solution in constant density regions
    when only oscillations (no interactions) are enabled.
)doc");

      class_object->def("Set_PositivityConstrain",&BaseType::Set_PositivityConstrain, py::arg("opt"),
R"doc(Toggle positivity constraint.

Parameters
----------
opt : bool
    If True, enforce positive flux during evolution.
)doc");

      class_object->def("Set_PositivityConstrainStep",&BaseType::Set_PositivityConstrainStep, py::arg("step"),
R"doc(Set positivity constraint step.

Parameters
----------
step : float
    Length scale for positivity enforcement.
)doc");

      class_object->def("Set_ProgressBar",&BaseType::Set_ProgressBar, py::arg("opt"),
R"doc(Toggle progress bar during evolution.

Parameters
----------
opt : bool
    If True, print progress bar to terminal.
)doc");

      // Mixing parameters
      class_object->def("Set_MixingParametersToDefault",&BaseType::Set_MixingParametersToDefault,
R"doc(Reset mixing parameters to default values.

Sets mixing angles, CP phases, and mass splittings to
values from Gonzalez-Garcia et al. global fits.
)doc");

      class_object->def("Set_Basis",&BaseType::Set_Basis, py::arg("basis"),
R"doc(Set the basis for evolution.

Parameters
----------
basis : Basis
    Evolution basis: mass or interaction (default).
    Interaction basis is recommended for most problems.
)doc");

      class_object->def("Set_MixingAngle",&BaseType::Set_MixingAngle,
        py::arg("i"),py::arg("j"),py::arg("angle"),
R"doc(Set a mixing angle.

Parameters
----------
i : int
    First state index (0-based).
j : int
    Second state index (must be > i).
angle : float
    Mixing angle in radians.

Example
-------
>>> nus.Set_MixingAngle(0, 1, 0.5904)  # theta_12
>>> nus.Set_MixingAngle(0, 2, 0.1496)  # theta_13
>>> nus.Set_MixingAngle(1, 2, 0.8553)  # theta_23
)doc");

      class_object->def("Get_MixingAngle",&BaseType::Get_MixingAngle,
        py::arg("i"),py::arg("j"),
R"doc(Get a mixing angle.

Parameters
----------
i : int
    First state index.
j : int
    Second state index.

Returns
-------
float
    Mixing angle in radians.
)doc");

      class_object->def("Set_CPPhase",&BaseType::Set_CPPhase,
        py::arg("i"),py::arg("j"),py::arg("phase"),
R"doc(Set a CP phase.

Parameters
----------
i : int
    First state index.
j : int
    Second state index.
phase : float
    CP phase in radians.

Example
-------
>>> nus.Set_CPPhase(0, 2, -1.38)  # delta_CP (standard 3-flavor)
)doc");

      class_object->def("Get_CPPhase",&BaseType::Get_CPPhase,
        py::arg("i"),py::arg("j"),
R"doc(Get a CP phase.

Parameters
----------
i : int
    First state index.
j : int
    Second state index.

Returns
-------
float
    CP phase in radians.
)doc");

      class_object->def("Set_SquareMassDifference",&BaseType::Set_SquareMassDifference,
        py::arg("i"),py::arg("dm2"),
R"doc(Set a squared mass difference.

Parameters
----------
i : int
    Mass eigenstate index (1, 2, ...).
dm2 : float
    Squared mass difference dm^2_{i1} in eV^2.

Example
-------
>>> nus.Set_SquareMassDifference(1, 7.42e-5)  # dm^2_21
>>> nus.Set_SquareMassDifference(2, 2.51e-3)  # dm^2_31
)doc");

      class_object->def("Get_SquareMassDifference",&BaseType::Get_SquareMassDifference,
        py::arg("i"),
R"doc(Get a squared mass difference.

Parameters
----------
i : int
    Mass eigenstate index.

Returns
-------
float
    Squared mass difference dm^2_{i1} in eV^2.
)doc");

      class_object->def("GetERange",&BaseType::GetERange,
R"doc(Get energy node values.

Returns
-------
array
    1D array of energy values in eV.
)doc");

      class_object->def("GetTrack",&BaseType::GetTrack,
R"doc(Get the current track object.

Returns
-------
Track
    The trajectory through the body.
)doc");

      class_object->def("GetBody",&BaseType::GetBody,
R"doc(Get the current body object.

Returns
-------
Body
    The propagation medium.
)doc");

      class_object->def("GetNumE",&BaseType::GetNumE,
R"doc(Get number of energy nodes.

Returns
-------
int
    Number of energy nodes.
)doc");

      class_object->def("GetNumRho",&BaseType::GetNumRho,
R"doc(Get number of density matrix equations.

Returns
-------
int
    1 if NT=neutrino or antineutrino, 2 if NT=both.
)doc");

      class_object->def("GetUseInteractions",&BaseType::GetUseInteractions,
R"doc(Check if interactions are enabled.

Returns
-------
bool
    True if non-coherent interactions are included.
)doc");

      class_object->def("GetUseOscillations",&BaseType::GetUseOscillations,
R"doc(Check if oscillations are enabled.

Returns
-------
bool
    True if oscillations are included.
)doc");

      class_object->def("InitializeInteractions",&BaseType::InitializeInteractions,
R"doc(Initialize interaction arrays.

Call this after setting up cross sections but before evolution
if you need to inspect the interaction structure.
Usually called automatically.
)doc");

      class_object->def("GetInteractionStructure",(std::shared_ptr<nusquids::nuSQUIDS::InteractionStructure>(BaseType::*)())&BaseType::GetInteractionStructure,
R"doc(Get the interaction structure.

Returns
-------
InteractionStructure
    Struct containing cross section tables and decay spectra.
    Can be reused for multiple nuSQUIDS instances.
)doc");

      class_object->def("GetHamiltonian",&BaseType::GetHamiltonian);

      class_object->def("GetTransformationMatrix",&BaseType::GetTransformationMatrix,
R"doc(Get the PMNS mixing matrix.

Returns
-------
gsl_matrix_complex
    The lepton mixing matrix U.
)doc");

      class_object->def("GetNeutrinoCrossSections",&BaseType::GetNeutrinoCrossSections,
R"doc(Get the neutrino cross section library.

Returns
-------
CrossSectionLibrary
    The cross section object used for interactions.
)doc");

      class_object->def("SetNeutrinoCrossSections",&BaseType::SetNeutrinoCrossSections,
        py::arg("xs"),
R"doc(Set a custom cross section library.

Parameters
----------
xs : CrossSectionLibrary
    Custom cross section library.
)doc");

      class_object->def("Set_Debug",&BaseType::Set_Debug, py::arg("debug"),
R"doc(Toggle debug output.

Parameters
----------
debug : bool
    If True, print debug information.
)doc");

      class_object->def("Set_IncludeOscillations",&BaseType::Set_IncludeOscillations);
      class_object->def("Set_GlashowResonance",&BaseType::Set_GlashowResonance);

      class_object->def("Set_NeutrinoSources",&BaseType::Set_NeutrinoSources, py::arg("opt"),
R"doc(Enable neutrino sources from bodies.

Parameters
----------
opt : bool
    If True, bodies can emit neutrinos.
)doc");

      class_object->def("Get_NeutrinoSources",&BaseType::Get_NeutrinoSources,
R"doc(Check if neutrino sources are enabled.

Returns
-------
bool
    True if body neutrino sources are enabled.
)doc");
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
      class_object = std::make_shared<py::class_<nuSQUIDSAtm<BaseType>, std::shared_ptr<nuSQUIDSAtm<BaseType>>>>(m,class_label.c_str(),
R"doc(Atmospheric neutrino propagation class.

nuSQUIDSAtm handles neutrino propagation through Earth for multiple
zenith angles simultaneously. It manages a collection of nuSQUIDS
objects, one per zenith angle, and provides methods to evaluate
fluxes as a function of both energy and cos(zenith).

This class is optimized for atmospheric neutrino experiments like
IceCube, ANTARES, Super-Kamiokande, etc.

Examples
--------
>>> import nuSQuIDS as nsq
>>> import numpy as np
>>> units = nsq.Const()
>>>
>>> # Set up energy and zenith grids
>>> E_nodes = nsq.logspace(100*units.GeV, 1*units.PeV, 100)
>>> cth_nodes = nsq.linspace(-1, 0.2, 20)
>>>
>>> # Create atmospheric propagator with interactions
>>> atm = nsq.nuSQUIDSAtm(cth_nodes, E_nodes, 3, nsq.NeutrinoType.both, True)
>>>
>>> # Set initial flux (shape: n_cth x n_E x 2 x n_flavors)
>>> # ... set flux array ...
>>> atm.Set_initial_state(flux, nsq.Basis.flavor)
>>> atm.EvolveState()
>>>
>>> # Evaluate flux
>>> phi = atm.EvalFlavor(1, -0.5, 1e11, 0)  # nu_mu at cth=-0.5, 100 GeV
)doc");

      class_object->def(py::init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType>(),
        py::arg("costh_nodes"),py::arg("E_nodes"),py::arg("numneu"),py::arg("NT"),
R"doc(Atmospheric mode constructor.

Parameters
----------
costh_nodes : array_like
    Cosine of zenith angle nodes (typically -1 to 0 for upgoing).
E_nodes : array_like
    Energy nodes in eV.
numneu : int
    Number of neutrino flavors.
NT : NeutrinoType
    neutrino, antineutrino, or both.
)doc");

      class_object->def(py::init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType,bool>(),
        py::arg("costh_nodes"),py::arg("E_nodes"),py::arg("numneu"),py::arg("NT"),py::arg("iinteraction"),
R"doc(Atmospheric mode constructor with interactions.

Parameters
----------
costh_nodes : array_like
    Cosine of zenith angle nodes.
E_nodes : array_like
    Energy nodes in eV.
numneu : int
    Number of neutrino flavors.
NT : NeutrinoType
    neutrino, antineutrino, or both.
iinteraction : bool
    If True, include non-coherent interactions.
)doc");

      class_object->def(py::init<marray<double,1>,marray<double,1>,unsigned int,NeutrinoType,bool,std::shared_ptr<CrossSectionLibrary>>(),
        py::arg("costh_nodes"),py::arg("E_nodes"),py::arg("numneu"),py::arg("NT"),py::arg("iinteraction"),py::arg("ncs"),
R"doc(Atmospheric mode constructor with custom cross sections.

Parameters
----------
costh_nodes : array_like
    Cosine of zenith angle nodes.
E_nodes : array_like
    Energy nodes in eV.
numneu : int
    Number of neutrino flavors.
NT : NeutrinoType
    neutrino, antineutrino, or both.
iinteraction : bool
    If True, include interactions.
ncs : CrossSectionLibrary
    Custom cross section library.
)doc");

      class_object->def(py::init<std::string>(),py::arg("filename"),
R"doc(Construct from HDF5 file.

Parameters
----------
filename : str
    Path to HDF5 file with saved nuSQUIDSAtm state.
)doc");

      class_object->def("EvolveState",&nuSQUIDSAtm<BaseType>::EvolveState,
R"doc(Evolve all zenith trajectories.

Propagates neutrinos through Earth for all zenith angles.
A progress bar is shown if Set_ProgressBar(True) was called.
)doc");

      class_object->def("Set_TauRegeneration",&nuSQUIDSAtm<BaseType>::Set_TauRegeneration, py::arg("opt"),
R"doc(Toggle tau regeneration for all trajectories.

Parameters
----------
opt : bool
    If True, include tau regeneration.
)doc");

      class_object->def("EvalFlavor",(double(nuSQUIDSAtm<BaseType>::*)(unsigned int,double,double,unsigned int,bool) const)&nuSQUIDSAtm<BaseType>::EvalFlavor,
        py::arg("flavor"),py::arg("costh"),py::arg("E"),py::arg("rho") = 0,py::arg("randomize_height") = false,
R"doc(Evaluate flux at given zenith and energy.

Parameters
----------
flavor : int
    Flavor index (0=nu_e, 1=nu_mu, 2=nu_tau).
costh : float
    Cosine of zenith angle.
E : float
    Energy in eV.
rho : int, optional
    Equation index: 0=neutrino (default), 1=antineutrino.
randomize_height : bool, optional
    If True, randomize production height (for systematics).

Returns
-------
float
    Interpolated flux value.
)doc");

      class_object->def("EvalFlavor",(double (nuSQUIDSAtm<BaseType>::*)(unsigned int, double, double, unsigned int, double, std::vector<bool>) const)&nuSQUIDSAtm<BaseType>::EvalFlavor,
        py::arg("flavor"),py::arg("costh"),py::arg("E"),py::arg("rho"),py::arg("scale"),py::arg("avr"),
R"doc(Evaluate flux with oscillation averaging.

Parameters
----------
flavor : int
    Flavor index.
costh : float
    Cosine of zenith.
E : float
    Energy in eV.
rho : int
    Equation index.
scale : float
    Averaging scale.
avr : list of bool
    Output: which scales were averaged.

Returns
-------
float
    Averaged flux value.
)doc");

      class_object->def("Set_EvalThreads",&nuSQUIDSAtm<BaseType>::Set_EvalThreads, py::arg("nthreads"),
R"doc(Set number of threads for evaluation.

Parameters
----------
nthreads : int
    Number of threads for parallel flux evaluation.
)doc");

      class_object->def("Get_EvalThreads",&nuSQUIDSAtm<BaseType>::Get_EvalThreads,
R"doc(Get number of evaluation threads.

Returns
-------
int
    Number of threads.
)doc");

      class_object->def("Set_EarthModel",&nuSQUIDSAtm<BaseType>::Set_EarthModel, py::arg("model"),
R"doc(Set Earth density model.

Parameters
----------
model : str
    Path to Earth model file (PREM format).
)doc");

      class_object->def("WriteStateHDF5",&nuSQUIDSAtm<BaseType>::WriteStateHDF5, py::arg("filename"), py::arg("overwrite") = true,
R"doc(Write state to HDF5 file.

Parameters
----------
filename : str
    Output filename.
overwrite : bool, optional
    If True (default), overwrite existing file.
)doc");

      class_object->def("ReadStateHDF5",&nuSQUIDSAtm<BaseType>::ReadStateHDF5, py::arg("filename"),
R"doc(Read state from HDF5 file.

Parameters
----------
filename : str
    Input filename.
)doc");

      class_object->def("Set_MixingAngle",&nuSQUIDSAtm<BaseType>::Set_MixingAngle,
        py::arg("i"),py::arg("j"),py::arg("angle"),
R"doc(Set mixing angle for all trajectories.

Parameters
----------
i : int
    First state index.
j : int
    Second state index.
angle : float
    Mixing angle in radians.
)doc");

      class_object->def("Get_MixingAngle",&nuSQUIDSAtm<BaseType>::Get_MixingAngle);

      class_object->def("Set_CPPhase",&nuSQUIDSAtm<BaseType>::Set_CPPhase,
        py::arg("i"),py::arg("j"),py::arg("phase"),
R"doc(Set CP phase for all trajectories.

Parameters
----------
i : int
    First state index.
j : int
    Second state index.
phase : float
    CP phase in radians.
)doc");

      class_object->def("Get_CPPhase",&nuSQUIDSAtm<BaseType>::Get_CPPhase);

      class_object->def("Set_SquareMassDifference",&nuSQUIDSAtm<BaseType>::Set_SquareMassDifference,
        py::arg("i"),py::arg("dm2"),
R"doc(Set squared mass difference for all trajectories.

Parameters
----------
i : int
    Mass eigenstate index.
dm2 : float
    Squared mass difference in eV^2.
)doc");

      class_object->def("Get_SquareMassDifference",&nuSQUIDSAtm<BaseType>::Get_SquareMassDifference);

      class_object->def("Set_h",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_h, py::arg("h"),
R"doc(Set step size for all trajectories.

Parameters
----------
h : float
    Step size in natural units.
)doc");

      class_object->def("Set_h",(void(nuSQUIDSAtm<BaseType>::*)(double,unsigned int))&nuSQUIDSAtm<BaseType>::Set_h,
        py::arg("h"),py::arg("ith"),
R"doc(Set step size for specific trajectory.

Parameters
----------
h : float
    Step size.
ith : int
    Trajectory index.
)doc");

      class_object->def("Set_h_max",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_h_max, py::arg("h_max"));
      class_object->def("Set_h_max",(void(nuSQUIDSAtm<BaseType>::*)(double,unsigned int))&nuSQUIDSAtm<BaseType>::Set_h_max);
      class_object->def("Set_h_min",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_h_min, py::arg("h_min"));
      class_object->def("Set_h_min",(void(nuSQUIDSAtm<BaseType>::*)(double,unsigned int))&nuSQUIDSAtm<BaseType>::Set_h_min);

      class_object->def("Set_ProgressBar",&nuSQUIDSAtm<BaseType>::Set_ProgressBar, py::arg("opt"),
R"doc(Toggle progress bar.

Parameters
----------
opt : bool
    If True, show progress during evolution.
)doc");

      class_object->def("Set_MixingParametersToDefault",&nuSQUIDSAtm<BaseType>::Set_MixingParametersToDefault,
R"doc(Reset mixing parameters to default values for all trajectories.)doc");

      class_object->def("Set_GSL_step",wrap_nusqatm_Set_GSL_STEP<BaseType>, py::arg("step_type"),
R"doc(Set GSL stepping algorithm for all trajectories.

Parameters
----------
step_type : GSL_STEP_FUNCTIONS
    Stepping algorithm (e.g., GSL_STEP_RKF45).
)doc");

      class_object->def("Set_rel_error",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_rel_error, py::arg("rel_error"),
R"doc(Set relative error for all trajectories.

Parameters
----------
rel_error : float
    Relative error tolerance.
)doc");

      class_object->def("Set_rel_error",(void(nuSQUIDSAtm<BaseType>::*)(double, unsigned int))&nuSQUIDSAtm<BaseType>::Set_rel_error);
      class_object->def("Set_abs_error",(void(nuSQUIDSAtm<BaseType>::*)(double))&nuSQUIDSAtm<BaseType>::Set_abs_error, py::arg("abs_error"));
      class_object->def("Set_abs_error",(void(nuSQUIDSAtm<BaseType>::*)(double, unsigned int))&nuSQUIDSAtm<BaseType>::Set_abs_error);
      class_object->def("Set_EvolLowPassCutoff",&nuSQUIDSAtm<BaseType>::Set_EvolLowPassCutoff);
      class_object->def("Set_EvolLowPassScale",&nuSQUIDSAtm<BaseType>::Set_EvolLowPassScale);

      class_object->def("GetNumE",&nuSQUIDSAtm<BaseType>::GetNumE,
R"doc(Get number of energy nodes.

Returns
-------
int
    Number of energy nodes.
)doc");

      class_object->def("GetNumCos",&nuSQUIDSAtm<BaseType>::GetNumCos,
R"doc(Get number of zenith nodes.

Returns
-------
int
    Number of cos(zenith) nodes.
)doc");

      class_object->def("GetNumNeu",&nuSQUIDSAtm<BaseType>::GetNumNeu,
R"doc(Get number of neutrino flavors.

Returns
-------
int
    Number of flavors.
)doc");

      class_object->def("GetNumRho",&nuSQUIDSAtm<BaseType>::GetNumRho,
R"doc(Get number of density matrix equations.

Returns
-------
int
    1 or 2 depending on NT.
)doc");

      class_object->def("GetnuSQuIDS",(std::vector<BaseType>&(nuSQUIDSAtm<BaseType>::*)())&nuSQUIDSAtm<BaseType>::GetnuSQuIDS,
        py::return_value_policy::reference_internal,
R"doc(Get all nuSQUIDS objects.

Returns
-------
list of nuSQUIDS
    The underlying nuSQUIDS objects for each zenith.
)doc");

      class_object->def("GetnuSQuIDS",(BaseType&(nuSQUIDSAtm<BaseType>::*)(unsigned int))&nuSQUIDSAtm<BaseType>::GetnuSQuIDS,
        py::return_value_policy::reference_internal, py::arg("ith"),
R"doc(Get nuSQUIDS object for specific zenith.

Parameters
----------
ith : int
    Zenith index.

Returns
-------
nuSQUIDS
    The nuSQUIDS object for that zenith angle.
)doc");

      class_object->def("Set_initial_state",(void(nuSQUIDSAtm<BaseType>::*)(const marray<double,3>&, Basis))&nuSQUIDSAtm<BaseType>::Set_initial_state,
        py::arg("ini_flux"),py::arg("basis") = Basis::flavor,
R"doc(Set initial flux (neutrino or antineutrino only).

Parameters
----------
ini_flux : array_like
    3D array of shape (n_costh, n_E, n_flavors).
basis : Basis, optional
    Basis of flux (default: flavor).
)doc");

      class_object->def("Set_initial_state",(void(nuSQUIDSAtm<BaseType>::*)(const marray<double,4>&, Basis))&nuSQUIDSAtm<BaseType>::Set_initial_state,
        py::arg("ini_flux"),py::arg("basis") = Basis::flavor,
R"doc(Set initial flux (both neutrino and antineutrino).

Parameters
----------
ini_flux : array_like
    4D array of shape (n_costh, n_E, 2, n_flavors).
    Axis 2: 0=neutrino, 1=antineutrino.
basis : Basis, optional
    Basis of flux (default: flavor).
)doc");

      class_object->def("GetStates", (marray<double,2>(nuSQUIDSAtm<BaseType>::*)(unsigned int))&nuSQUIDSAtm<BaseType>::GetStates,
        py::arg("rho") = 0,
R"doc(Get evolved states of all nodes.

Parameters
----------
rho : int, optional
    Equation index (default: 0).

Returns
-------
array
    2D array of states.
)doc");

      class_object->def("GetERange",&nuSQUIDSAtm<BaseType>::GetERange,
R"doc(Get energy node values.

Returns
-------
array
    1D array of energies in eV.
)doc");

      class_object->def("GetCosthRange",&nuSQUIDSAtm<BaseType>::GetCosthRange,
R"doc(Get cos(zenith) node values.

Returns
-------
array
    1D array of cos(zenith) values.
)doc");

      class_object->def("Set_IncludeOscillations",&nuSQUIDSAtm<BaseType>::Set_IncludeOscillations, py::arg("opt"));
      class_object->def("Set_GlashowResonance",&nuSQUIDSAtm<BaseType>::Set_GlashowResonance, py::arg("opt"));
      class_object->def("Set_TauRegeneration",&nuSQUIDSAtm<BaseType>::Set_TauRegeneration);
      class_object->def("Set_AllowConstantDensityOscillationOnlyEvolution",&nuSQUIDSAtm<BaseType>::Set_AllowConstantDensityOscillationOnlyEvolution, py::arg("opt"));
      class_object->def("Set_PositivyConstrain",&nuSQUIDSAtm<BaseType>::Set_PositivityConstrain, py::arg("opt"));
      class_object->def("Set_PositivyConstrainStep",&nuSQUIDSAtm<BaseType>::Set_PositivityConstrainStep, py::arg("step"));
      class_object->def("Get_EvalThreads",&nuSQUIDSAtm<BaseType>::Get_EvalThreads);
      class_object->def("Set_EvalThreads",&nuSQUIDSAtm<BaseType>::Set_EvalThreads);
      class_object->def("Set_EarthModel",&nuSQUIDSAtm<BaseType>::Set_EarthModel);
      class_object->def("SetNeutrinoCrossSections",&nuSQUIDSAtm<BaseType>::SetNeutrinoCrossSections, py::arg("xs"));
      class_object->def("GetNeutrinoCrossSections",&nuSQUIDSAtm<BaseType>::GetNeutrinoCrossSections);
    }
    std::shared_ptr<py::class_<nuSQUIDSAtm<BaseType>, std::shared_ptr<nuSQUIDSAtm<BaseType>>>> GetClassObject() {
      return class_object;
    }
};

#endif
