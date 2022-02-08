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

#include "nuSQUIDSpy.h"

PYBIND11_MODULE(nuSQuIDS, m)
{
  m.doc() = "nuSQuIDS Python Bindings"; // module docstring
  // import numpy array definitions
  //np::initialize();
  import_ufunc();
#if PY_VERSION_HEX >= 0x03000000
  import_array1();
#else
  import_array();
#endif

  py::enum_<GSL_STEP_FUNCTIONS>(m,"GSL_STEP_FUNCTIONS")
    .value("GSL_STEP_RK2",GSL_STEP_RK2)
    .value("GSL_STEP_RK4",GSL_STEP_RK4)
    .value("GSL_STEP_RKF45",GSL_STEP_RKF45)
    .value("GSL_STEP_RKCK",GSL_STEP_RKCK)
    .value("GSL_STEP_RK8PD",GSL_STEP_RK8PD)
    .value("GSL_STEP_MSADAMS",GSL_STEP_MSADAMS)
    .export_values()
  ;

  py::enum_<Basis>(m,"Basis")
    .value("mass",mass)
    .value("flavor",flavor)
    .value("interaction",interaction)
    .export_values()
  ;

  py::enum_<NeutrinoCrossSections::NeutrinoFlavor>(m,"NeutrinoCrossSections_NeutrinoFlavor")
    .value("electron",NeutrinoCrossSections::NeutrinoFlavor::electron)
    .value("muon",NeutrinoCrossSections::NeutrinoFlavor::muon)
    .value("tau",NeutrinoCrossSections::NeutrinoFlavor::tau)
    .value("sterile",NeutrinoCrossSections::NeutrinoFlavor::sterile)
    .export_values()
  ;

  py::enum_<NeutrinoCrossSections::NeutrinoType>(m,"NeutrinoCrossSections_NeutrinoType")
    .value("neutrino",NeutrinoCrossSections::NeutrinoType::neutrino)
    .value("antineutrino",NeutrinoCrossSections::NeutrinoType::antineutrino)
    .export_values()
  ;

  py::enum_<NeutrinoCrossSections::Current>(m,"NeutrinoCrossSections_Current")
    .value("CC",NeutrinoCrossSections::Current::CC)
    .value("NC",NeutrinoCrossSections::Current::NC)
    .value("GR",NeutrinoCrossSections::Current::GR)
    .export_values()
  ;


  py::class_<squids::SU_vector, std::shared_ptr<squids::SU_vector>>(m,"SU_vector")
    .def(py::init< std::vector<double>>())
    .def(py::init<unsigned int>())
    .def(py::init<const squids::SU_vector&>())
    .def("Dim",&squids::SU_vector::Dim)
    .def("Size",&squids::SU_vector::Size)
    .def("SetAllComponents",&squids::SU_vector::SetAllComponents)
    .def("GetComponents",&squids::SU_vector::GetComponents)
    .def("Rotate",(squids::SU_vector(squids::SU_vector::*)(unsigned int, unsigned int, double, double) const)&squids::SU_vector::Rotate)
    .def("RotateToB0",&squids::SU_vector::RotateToB0)
    .def("RotateToB1",&squids::SU_vector::RotateToB1)
    .def("WeightedRotation",(void(squids::SU_vector::*)(const squids::Const&, const squids::SU_vector&, const squids::Const&))&squids::SU_vector::WeightedRotation)
    .def("Transpose",&squids::SU_vector::Transpose)
    .def("Imag",&squids::SU_vector::Imag)
    .def("Real",&squids::SU_vector::Real)
    //.def("UTransform",(squids::SU_vector(squids::SU_vector::*)(const squids::SU_vector, gsl_complex))&squids::SU_vector::UTransform)
    //.def("GetEigenSystem",(std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)>(squids::SU_vector::*)())&squids::SU_vector::GetEigenSystem)
    //.def("GetGSLMatrix",(std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)>(squids::SU_vector::*)())&squids::SU_vector::GetGSLMatrix)
    .def(py::self += py::self)
    .def(py::self + py::self)
    .def(py::self -= py::self)
    .def(py::self - py::self)
    .def(py::self * double())
    .def(double() * py::self)
    .def(py::self *= double())
    .def(py::self * py::self)
    .def(py::self /= double())
    .def(py::self == py::self)
    .def(-py::self)
    .def("__str__",PrintObject<squids::SU_vector>)
    .def("Projector",&squids::SU_vector::Projector)
    .def_static("Identity",&squids::SU_vector::Identity)
    .def("PosProjector",&squids::SU_vector::PosProjector)
    .def("NegProjector",&squids::SU_vector::NegProjector)
    .def_static("Generator",&squids::SU_vector::Generator)
  ;

  py::class_<squids::detail::AdditionProxy, std::shared_ptr<squids::detail::AdditionProxy>>(m,"AdditionProxy");
  py::implicitly_convertible<squids::detail::AdditionProxy , squids::SU_vector>();

  py::class_<squids::detail::SubtractionProxy, std::shared_ptr<squids::detail::SubtractionProxy>>(m,"SubtractionProxy");
  py::implicitly_convertible<squids::detail::SubtractionProxy, squids::SU_vector>();

  py::class_<squids::detail::NegationProxy, std::shared_ptr<squids::detail::NegationProxy>>(m,"NegationProxy");
  py::implicitly_convertible<squids::detail::NegationProxy, squids::SU_vector>();

  py::class_<squids::detail::EvolutionProxy, std::shared_ptr<squids::detail::EvolutionProxy>>(m,"EvolutionProxy");
  py::implicitly_convertible<squids::detail::EvolutionProxy, squids::SU_vector>();

  py::class_<squids::detail::MultiplicationProxy, std::shared_ptr<squids::detail::MultiplicationProxy>>(m,"MultiplicationProxy");
  py::implicitly_convertible<squids::detail::MultiplicationProxy, squids::SU_vector>();

  py::class_<squids::detail::iCommutatorProxy, std::shared_ptr<squids::detail::iCommutatorProxy>>(m,"iCommutatorProxy");
  py::implicitly_convertible<squids::detail::iCommutatorProxy, squids::SU_vector>();

  py::class_<squids::detail::ACommutatorProxy, std::shared_ptr<squids::detail::ACommutatorProxy>>(m,"ACommutatorProxy");
  py::implicitly_convertible<squids::detail::ACommutatorProxy, squids::SU_vector>();

  py::enum_<NeutrinoType>(m,"NeutrinoType")
    .value("neutrino",neutrino)
    .value("antineutrino",antineutrino)
    .value("both",both)
    .export_values()
  ;

  m.def("linspace",linspace,py::arg("min"),py::arg("max"),py::arg("samples"));
  m.def("logspace",logspace,py::arg("min"),py::arg("max"),py::arg("samples"));

  RegisterBasicNuSQuIDSPythonBindings<nuSQUIDS>(m,"nuSQUIDS");
  RegisterBasicAtmNuSQuIDSPythonBindings<nuSQUIDS>(m,"nuSQUIDSAtm");

  py::class_<NeutrinoCrossSections, std::shared_ptr<NeutrinoCrossSections>>(m,"NeutrinoCrossSections");

  py::class_<NullCrossSections, NeutrinoCrossSections, std::shared_ptr<NullCrossSections>>(m,"NullCrossSections")
    .def("TotalCrossSection",&NullCrossSections::TotalCrossSection)
    .def("SingleDifferentialCrossSection",&NullCrossSections::SingleDifferentialCrossSection)
  ;

  py::class_<GlashowResonanceCrossSection, NeutrinoCrossSections, std::shared_ptr<GlashowResonanceCrossSection>>(m,"GlashowResonanceCrossSection")
    .def("TotalCrossSection",&GlashowResonanceCrossSection::TotalCrossSection)
    .def("SingleDifferentialCrossSection",&GlashowResonanceCrossSection::SingleDifferentialCrossSection)
    .def("WDecayBranchingFraction",&GlashowResonanceCrossSection::WDecayBranchingFraction)
  ;

  py::class_<NeutrinoDISCrossSectionsFromTables, NeutrinoCrossSections, std::shared_ptr<NeutrinoDISCrossSectionsFromTables>>(m,"NeutrinoDISCrossSectionsFromTables")
    .def(py::init<std::string>())
    .def("TotalCrossSection",&NeutrinoDISCrossSectionsFromTables::TotalCrossSection)
    .def("SingleDifferentialCrossSection",&NeutrinoDISCrossSectionsFromTables::SingleDifferentialCrossSection)
    .def("GetEmin",&NeutrinoDISCrossSectionsFromTables::GetEmin)
    .def("GetEmax",&NeutrinoDISCrossSectionsFromTables::GetEmax)
    .def("WriteHDF",&NeutrinoDISCrossSectionsFromTables::WriteHDF)
    .def("WriteText",&NeutrinoDISCrossSectionsFromTables::WriteText)
  ;

  py::enum_<PDGCode>(m,"PDGCode")
    .value("electron",electron)
    .value("isoscalar_nucleon",isoscalar_nucleon)
    .value("proton",proton)
    .value("neutron",neutron)
    .export_values()
  ;

  py::class_<CrossSectionLibrary, std::shared_ptr<CrossSectionLibrary>>(m,"CrossSectionLibrary")
    .def(py::init<>())
    //TODO: map constructor?
    .def("crossSectionForTarget", &CrossSectionLibrary::crossSectionForTarget)
    .def("hasTarget", (bool(CrossSectionLibrary::*)(typename std::underlying_type<PDGCode>::type))&CrossSectionLibrary::hasTarget)
    .def("addTarget", (void(CrossSectionLibrary::*)(typename std::underlying_type<PDGCode>::type, std::shared_ptr<NeutrinoCrossSections>))&CrossSectionLibrary::hasTarget)
  ;

  // what is this?
  //bp::def("loadDefaultCrossSections",loadDefaultCrossSections);

  py::class_<TauDecaySpectra, std::shared_ptr<TauDecaySpectra>>(m,"TauDecaySpectra")
    .def(py::init<marray<double,1>>())
    .def("dNdEnu_All",&TauDecaySpectra::dNdEnu_All)
    .def("dNdEnu_Lep",&TauDecaySpectra::dNdEnu_Lep)
    .def("dNdEle_All",&TauDecaySpectra::dNdEle_All)
    .def("dNdEle_Lep",&TauDecaySpectra::dNdEle_Lep)
    .def("GetTauToHadronBranchingRatio",&TauDecaySpectra::GetTauToHadronBranchingRatio)
    .def("GetTauToLeptonBranchingRatio",&TauDecaySpectra::GetTauToLeptonBranchingRatio)
    .def("TauDecayToLepton",&TauDecaySpectra::TauDecayToLepton)
    .def("TauDecayToHadron",&TauDecaySpectra::TauDecayToHadron)
    .def("TauDecayToAll",&TauDecaySpectra::TauDecayToAll)
    .def("TauDecayToPion",&TauDecaySpectra::TauDecayToPion)
    .def("TauDecayToRho",&TauDecaySpectra::TauDecayToRho)
    .def("TauDecayToA1",&TauDecaySpectra::TauDecayToA1)
  ;

  py::class_<squids::Const>(m,"Const")
    .def_readonly("GF",&squids::Const::GF)
    .def_readonly("Na",&squids::Const::Na)
    .def_readonly("GF",&squids::Const::GF)
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
    auto outer
    = py::class_<Body, std::shared_ptr<Body>>(m,"Body")
    .def("density",&Body::density)
    .def("ye",&Body::ye)
    ;

    py::class_<Body::Track, std::shared_ptr<Body::Track>>(outer,"Track")
    .def("GetInitialX",&Body::Track::GetInitialX)
    .def("GetFinalX",&Body::Track::GetFinalX)
    .def("GetX",&Body::Track::GetX)
    .def("SetX",&Body::Track::SetX)
    .def("ReverseTrack",&Body::Track::ReverseTrack)
    ;
  }

  {
    auto outer
    = py::class_<Vacuum, Body, std::shared_ptr<Vacuum>>(m,"Vacuum")
    .def(py::init<>())
    .def("density",&Vacuum::density,py::arg("track"))
    .def("ye",&Vacuum::ye,py::arg("track"))
    ;

    py::class_<Vacuum::Track, Body::Track, std::shared_ptr<Vacuum::Track>>(outer,"Track")
    .def(py::init<double>(),py::arg("xend"))
    .def(py::init<double,double>(),py::arg("xini"),py::arg("xend"))
    .def(py::init<double,double,double>(),py::arg("x"),py::arg("xini"),py::arg("xend"))
    .def("GetInitialX",&Vacuum::Track::GetInitialX)
    .def("GetFinalX",&Vacuum::Track::GetFinalX)
    .def("GetX",&Vacuum::Track::GetX)
    .def("SetX",&Vacuum::Track::SetX)
    .def("ReverseTrack",&Vacuum::Track::ReverseTrack)
    ;
  }

  {
    auto outer
    = py::class_<ConstantDensity, Body, std::shared_ptr<ConstantDensity>>(m,"ConstantDensity")
    .def(py::init<double,double>())
    .def("density",&ConstantDensity::density)
    .def("ye",&ConstantDensity::ye)
    ;

    py::class_<ConstantDensity::Track, Body::Track, std::shared_ptr<ConstantDensity::Track>>(outer,"Track")
    .def(py::init<double>())
    .def(py::init<double,double>())
    .def("GetInitialX",&ConstantDensity::Track::GetInitialX)
    .def("GetFinalX",&ConstantDensity::Track::GetFinalX)
    .def("GetX",&ConstantDensity::Track::GetX)
    .def("SetX",&ConstantDensity::Track::SetX)
    .def("ReverseTrack",&ConstantDensity::Track::ReverseTrack)
    ;
  }

  {
    auto outer
    = py::class_<VariableDensity, Body, std::shared_ptr<VariableDensity>>(m,"VariableDensity")
    .def(py::init< std::vector<double>,std::vector<double>,std::vector<double>>())
    ;

    py::class_<VariableDensity::Track, Body::Track, std::shared_ptr<VariableDensity::Track>>(outer,"Track")
    .def(py::init<double>())
    .def(py::init<double,double>())
    .def("GetInitialX",&VariableDensity::Track::GetInitialX)
    .def("GetFinalX",&VariableDensity::Track::GetFinalX)
    .def("GetX",&VariableDensity::Track::GetX)
    .def("SetX",&VariableDensity::Track::SetX)
    .def("ReverseTrack",&VariableDensity::Track::ReverseTrack)
    ;
  }

  {
    auto outer
    = py::class_<Earth, Body, std::shared_ptr<Earth>>(m,"Earth")
    .def(py::init<std::string>())
    ;

    py::class_<Earth::Track, Body::Track, std::shared_ptr<Earth::Track>>(outer,"Track")
    .def(py::init<double>())
    .def(py::init<double,double,double>())
    .def("GetInitialX",&Earth::Track::GetInitialX)
    .def("GetFinalX",&Earth::Track::GetFinalX)
    .def("GetX",&Earth::Track::GetX)
    .def("SetX",&Earth::Track::SetX)
    .def("ReverseTrack",&Earth::Track::ReverseTrack)
    ;
  }

  {
    auto outer
    = py::class_<Sun, Body, std::shared_ptr<Sun>>(m,"Sun")
    .def(py::init<std::string>())
    ;

    py::class_<Sun::Track, Body::Track, std::shared_ptr<Sun::Track>>(outer,"Track")
    .def(py::init<double>())
    .def(py::init<double,double>())
    .def("GetInitialX",&Sun::Track::GetInitialX)
    .def("GetFinalX",&Sun::Track::GetFinalX)
    .def("GetX",&Sun::Track::GetX)
    .def("SetX",&Sun::Track::SetX)
    .def("ReverseTrack",&Sun::Track::ReverseTrack)
    ;
  }

  {
    auto outer
    = py::class_<SunASnu, Body, std::shared_ptr<SunASnu>>(m,"SunASnu")
    .def(py::init<std::string>())
    ;

    py::class_<SunASnu::Track, Body::Track, std::shared_ptr<SunASnu::Track>>(outer,"Track")
    .def(py::init<double>())
    .def(py::init<double,double>())
    .def("GetInitialX",&SunASnu::Track::GetInitialX)
    .def("GetFinalX",&SunASnu::Track::GetFinalX)
    .def("GetX",&SunASnu::Track::GetX)
    .def("SetX",&SunASnu::Track::SetX)
    .def("ReverseTrack",&SunASnu::Track::ReverseTrack)
    ;
  }

  {
    auto outer
    = py::class_<EarthAtm, Body, std::shared_ptr<EarthAtm>>(m,"EarthAtm")
    .def(py::init<std::string>())
    .def("GetRadius",&EarthAtm::GetRadius)
    .def("GetAtmosphereHeight",&EarthAtm::GetAtmosphereHeight)
    .def("SetAtmosphereHeight",&EarthAtm::GetAtmosphereHeight)
    .def("MakeTrack",&EarthAtm::MakeTrack)
    .def("MakeTrackWithCosine",&EarthAtm::MakeTrackWithCosine)
    ;

    py::class_<EarthAtm::Track, Body::Track, std::shared_ptr<EarthAtm::Track>>(outer,"Track")
    .def(py::init<EarthAtm::Track>())
    .def(py::init<double,double,double>())
    .def(py::init<double,double,double,double>())
    .def("GetInitialX",&EarthAtm::Track::GetInitialX)
    .def("GetFinalX",&EarthAtm::Track::GetFinalX)
    .def("GetX",&EarthAtm::Track::GetX)
    .def("SetX",&EarthAtm::Track::SetX)
    .def("ReverseTrack",&EarthAtm::Track::ReverseTrack)
    ;
  }

  /*
  from_python_sequence<std::vector<double>, variable_capacity_policy>();
  to_python_converter<std::vector<double, class std::allocator<double>>, VecToList<double>> ();
  // register marray converters
  py::to_python_converter<marray<double,1> , marray_to_numpyarray<1>>();
  py::to_python_converter<marray<double,2> , marray_to_numpyarray<2>>();
  py::to_python_converter<marray<double,3> , marray_to_numpyarray<3>>();
  py::to_python_converter<marray<double,4> , marray_to_numpyarray<4>>();

  marray_from_python<1>();
  marray_from_python<2>();
  marray_from_python<3>();
  marray_from_python<4>();
  */
}
