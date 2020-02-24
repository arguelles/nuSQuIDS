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

BOOST_PYTHON_MODULE(nuSQUIDSpy)
{
  // import numpy array definitions
  //np::initialize();
  import_ufunc();
#if PY_VERSION_HEX >= 0x03000000
  import_array1();
#else
  import_array();
#endif

  enum_<GSL_STEP_FUNCTIONS>("GSL_STEP_FUNCTIONS")
    .value("GSL_STEP_RK2",GSL_STEP_RK2)
    .value("GSL_STEP_RK4",GSL_STEP_RK4)
    .value("GSL_STEP_RKF45",GSL_STEP_RKF45)
    .value("GSL_STEP_RKCK",GSL_STEP_RKCK)
    .value("GSL_STEP_RK8PD",GSL_STEP_RK8PD)
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

  class_<squids::SU_vector, std::shared_ptr<squids::SU_vector> >("SU_vector")
    .def(init< std::vector<double> >())
    .def(init<unsigned int>())
    .def(init<const squids::SU_vector&>())
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
    .def(self += self)
    .def(self + self)
    .def(self -= self)
    .def(self - self)
    .def(self * double())
    .def(double() * self)
    .def(self *= double())
    .def(self * self)
    .def(self /= double())
    .def(self == self)
    .def(-self)
    .def(self_ns::str(self_ns::self))
    .def("Projector",&squids::SU_vector::Projector)
    .def("Identity",&squids::SU_vector::Identity).staticmethod("Identity")
    .def("PosProjector",&squids::SU_vector::PosProjector)
    .def("NegProjector",&squids::SU_vector::NegProjector)
    .def("Generator",&squids::SU_vector::Generator).staticmethod("Generator")
  ;

  class_<squids::detail::AdditionProxy, std::shared_ptr<squids::detail::AdditionProxy>>("AdditionProxy", no_init);
  implicitly_convertible< squids::detail::AdditionProxy , squids::SU_vector >();

  class_<squids::detail::SubtractionProxy, std::shared_ptr<squids::detail::SubtractionProxy>>("SubtractionProxy", no_init);
  implicitly_convertible< squids::detail::SubtractionProxy, squids::SU_vector >();

  class_<squids::detail::NegationProxy, std::shared_ptr<squids::detail::NegationProxy>>("NegationProxy", no_init);
  implicitly_convertible< squids::detail::NegationProxy, squids::SU_vector >();

  class_<squids::detail::EvolutionProxy, std::shared_ptr<squids::detail::EvolutionProxy>>("EvolutionProxy", no_init);
  implicitly_convertible< squids::detail::EvolutionProxy, squids::SU_vector >();

  class_<squids::detail::MultiplicationProxy, std::shared_ptr<squids::detail::MultiplicationProxy>>("MultiplicationProxy", no_init);
  implicitly_convertible< squids::detail::MultiplicationProxy, squids::SU_vector >();

  class_<squids::detail::iCommutatorProxy, std::shared_ptr<squids::detail::iCommutatorProxy>>("iCommutatorProxy", no_init);
  implicitly_convertible< squids::detail::iCommutatorProxy, squids::SU_vector >();

  class_<squids::detail::ACommutatorProxy, std::shared_ptr<squids::detail::ACommutatorProxy>>("ACommutatorProxy", no_init);
  implicitly_convertible< squids::detail::ACommutatorProxy, squids::SU_vector >();

  enum_<NeutrinoType>("NeutrinoType")
    .value("neutrino",neutrino)
    .value("antineutrino",antineutrino)
    .value("both",both)
  ;

  bp::def("linspace",linspace,bp::args("min","max","samples"));
  bp::def("logspace",logspace,bp::args("min","max","samples"));

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

  class_<TauDecaySpectra, std::shared_ptr<TauDecaySpectra>, boost::noncopyable>("TauDecaySpectra")
    .def(init<marray<double,1>>())
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
