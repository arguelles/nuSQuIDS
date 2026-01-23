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

#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/resources.h>
#include <nuSQuIDS/marray.h>

#include <string>
#include <vector>
#include <memory>
#include <sstream>

using namespace nusquids;

// Helper to convert std::vector to marray<double,1>
marray<double, 1> vec_to_marray1(const std::vector<double>& vec) {
    marray<double, 1> arr;
    arr.resize(std::array<size_t,1>{vec.size()});
    std::copy(vec.begin(), vec.end(), arr.get_data());
    return arr;
}

// Helper to convert 2D Julia array (as flat vector + dims) to marray<double,2>
marray<double, 2> vec_to_marray2(const std::vector<double>& vec, size_t dim1, size_t dim2) {
    marray<double, 2> arr;
    arr.resize(std::array<size_t,2>{dim1, dim2});
    std::copy(vec.begin(), vec.end(), arr.get_data());
    return arr;
}

// Helper to convert 3D Julia array (as flat vector + dims) to marray<double,3>
marray<double, 3> vec_to_marray3(const std::vector<double>& vec, size_t dim1, size_t dim2, size_t dim3) {
    marray<double, 3> arr;
    arr.resize(std::array<size_t,3>{dim1, dim2, dim3});
    std::copy(vec.begin(), vec.end(), arr.get_data());
    return arr;
}

// Helper function to create energy arrays
std::vector<double> linspace_vec(double min, double max, size_t samples) {
    std::vector<double> result(samples);
    double step = (max - min) / (samples - 1);
    for (size_t i = 0; i < samples; ++i) {
        result[i] = min + i * step;
    }
    return result;
}

std::vector<double> logspace_vec(double min, double max, size_t samples) {
    std::vector<double> result(samples);
    double log_min = std::log10(min);
    double log_max = std::log10(max);
    double step = (log_max - log_min) / (samples - 1);
    for (size_t i = 0; i < samples; ++i) {
        result[i] = std::pow(10.0, log_min + i * step);
    }
    return result;
}

// GSL step function enum (mirrors Python bindings)
enum GSL_STEP_FUNCTIONS {
    GSL_STEP_RK2,
    GSL_STEP_RK4,
    GSL_STEP_RKF45,
    GSL_STEP_RKCK,
    GSL_STEP_RK8PD,
    GSL_STEP_MSADAMS
};

template<typename BaseType>
void wrap_Set_GSL_STEP(BaseType* nusq, GSL_STEP_FUNCTIONS step_enum) {
    switch(step_enum) {
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

// Declare supertype relationships for inheritance
namespace jlcxx {
    // Body hierarchy
    template<> struct SuperType<Vacuum> { typedef Body type; };
    template<> struct SuperType<ConstantDensity> { typedef Body type; };
    template<> struct SuperType<VariableDensity> { typedef Body type; };
    template<> struct SuperType<Earth> { typedef Body type; };
    template<> struct SuperType<EarthAtm> { typedef Body type; };
    template<> struct SuperType<Sun> { typedef Body type; };
    template<> struct SuperType<SunASnu> { typedef Body type; };

    // Track hierarchy
    template<> struct SuperType<Vacuum::Track> { typedef Body::Track type; };
    template<> struct SuperType<ConstantDensity::Track> { typedef Body::Track type; };
    template<> struct SuperType<VariableDensity::Track> { typedef Body::Track type; };
    template<> struct SuperType<Earth::Track> { typedef Body::Track type; };
    template<> struct SuperType<EarthAtm::Track> { typedef Body::Track type; };
    template<> struct SuperType<Sun::Track> { typedef Body::Track type; };
    template<> struct SuperType<SunASnu::Track> { typedef Body::Track type; };

    // Cross section hierarchy
    template<> struct SuperType<NullCrossSections> { typedef NeutrinoCrossSections type; };
    template<> struct SuperType<GlashowResonanceCrossSection> { typedef NeutrinoCrossSections type; };
    template<> struct SuperType<NeutrinoDISCrossSectionsFromTables> { typedef NeutrinoCrossSections type; };
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    // ========== Enumerations ==========

    mod.add_bits<Basis>("Basis", jlcxx::julia_type("CppEnum"));
    mod.set_const("mass", mass);
    mod.set_const("flavor", flavor);
    mod.set_const("interaction", interaction);

    mod.add_bits<NeutrinoType>("NeutrinoType", jlcxx::julia_type("CppEnum"));
    mod.set_const("neutrino", neutrino);
    mod.set_const("antineutrino", antineutrino);
    mod.set_const("both", both);

    mod.add_bits<GSL_STEP_FUNCTIONS>("GSLStepFunction", jlcxx::julia_type("CppEnum"));
    mod.set_const("GSL_STEP_RK2", GSL_STEP_RK2);
    mod.set_const("GSL_STEP_RK4", GSL_STEP_RK4);
    mod.set_const("GSL_STEP_RKF45", GSL_STEP_RKF45);
    mod.set_const("GSL_STEP_RKCK", GSL_STEP_RKCK);
    mod.set_const("GSL_STEP_RK8PD", GSL_STEP_RK8PD);
    mod.set_const("GSL_STEP_MSADAMS", GSL_STEP_MSADAMS);

    mod.add_bits<NeutrinoCrossSections::NeutrinoFlavor>("NeutrinoFlavor", jlcxx::julia_type("CppEnum"));
    mod.set_const("electron_flavor", NeutrinoCrossSections::NeutrinoFlavor::electron);
    mod.set_const("muon_flavor", NeutrinoCrossSections::NeutrinoFlavor::muon);
    mod.set_const("tau_flavor", NeutrinoCrossSections::NeutrinoFlavor::tau);
    mod.set_const("sterile_flavor", NeutrinoCrossSections::NeutrinoFlavor::sterile);

    mod.add_bits<NeutrinoCrossSections::Current>("Current", jlcxx::julia_type("CppEnum"));
    mod.set_const("CC", NeutrinoCrossSections::Current::CC);
    mod.set_const("NC", NeutrinoCrossSections::Current::NC);
    mod.set_const("GR", NeutrinoCrossSections::Current::GR);

    mod.add_bits<PDGCode>("PDGCode", jlcxx::julia_type("CppEnum"));
    mod.set_const("electron_pdg", electron);
    mod.set_const("isoscalar_nucleon_pdg", isoscalar_nucleon);
    mod.set_const("proton_pdg", proton);
    mod.set_const("neutron_pdg", neutron);
    mod.set_const("hydrogen_pdg", hydrogen);
    mod.set_const("deuteron_pdg", deuteron);
    mod.set_const("helium3_pdg", helium3);
    mod.set_const("helium4_pdg", helium4);
    mod.set_const("carbon_pdg", carbon);
    mod.set_const("nitrogen_pdg", nitrogen);
    mod.set_const("oxygen_pdg", oxygen);
    mod.set_const("sodium_pdg", sodium);
    mod.set_const("magnesium_pdg", magnesium);
    mod.set_const("aluminum_pdg", aluminum);
    mod.set_const("silicon_pdg", silicon);
    mod.set_const("sulfur_pdg", sulfur);
    mod.set_const("calcium_pdg", calcium);
    mod.set_const("iron_pdg", iron);
    mod.set_const("nickel_pdg", nickel);
    mod.set_const("tungsten_pdg", tungsten);
    mod.set_const("lead_pdg", lead);

    // ========== Physical Constants (squids::Const) ==========

    mod.add_type<squids::Const>("Const")
        .constructor<>()
        .method("GF", [](const squids::Const& c) { return c.GF; })
        .method("Na", [](const squids::Const& c) { return c.Na; })
        .method("PeV", [](const squids::Const& c) { return c.PeV; })
        .method("TeV", [](const squids::Const& c) { return c.TeV; })
        .method("GeV", [](const squids::Const& c) { return c.GeV; })
        .method("MeV", [](const squids::Const& c) { return c.MeV; })
        .method("keV", [](const squids::Const& c) { return c.keV; })
        .method("eV", [](const squids::Const& c) { return c.eV; })
        .method("kg", [](const squids::Const& c) { return c.kg; })
        .method("gr", [](const squids::Const& c) { return c.gr; })
        .method("meter", [](const squids::Const& c) { return c.meter; })
        .method("cm", [](const squids::Const& c) { return c.cm; })
        .method("km", [](const squids::Const& c) { return c.km; })
        .method("fermi", [](const squids::Const& c) { return c.fermi; })
        .method("angstrom", [](const squids::Const& c) { return c.angstrom; })
        .method("AU", [](const squids::Const& c) { return c.AU; })
        .method("parsec", [](const squids::Const& c) { return c.parsec; })
        .method("picobarn", [](const squids::Const& c) { return c.picobarn; })
        .method("femtobarn", [](const squids::Const& c) { return c.femtobarn; })
        .method("sec", [](const squids::Const& c) { return c.sec; })
        .method("hour", [](const squids::Const& c) { return c.hour; })
        .method("day", [](const squids::Const& c) { return c.day; })
        .method("year", [](const squids::Const& c) { return c.year; });

    // ========== SU_vector ==========

    mod.add_type<squids::SU_vector>("SU_vector")
        .constructor<unsigned int>()
        .constructor<std::vector<double>>()
        .method("Dim", &squids::SU_vector::Dim)
        .method("Size", &squids::SU_vector::Size)
        .method("SetAllComponents", &squids::SU_vector::SetAllComponents)
        .method("GetComponents", &squids::SU_vector::GetComponents)
        .method("Transpose", &squids::SU_vector::Transpose)
        .method("Imag", &squids::SU_vector::Imag)
        .method("Real", &squids::SU_vector::Real);

    // Static methods for SU_vector (take dimension as first argument)
    mod.method("SU_vector_Projector", &squids::SU_vector::Projector);
    mod.method("SU_vector_PosProjector", &squids::SU_vector::PosProjector);
    mod.method("SU_vector_NegProjector", &squids::SU_vector::NegProjector);
    mod.method("SU_vector_Identity", &squids::SU_vector::Identity);
    mod.method("SU_vector_Generator", &squids::SU_vector::Generator);

    // ========== Body and Track Base Classes ==========

    mod.add_type<Body::Track>("BodyTrack")
        .method("GetInitialX", &Body::Track::GetInitialX)
        .method("GetFinalX", &Body::Track::GetFinalX)
        .method("GetX", &Body::Track::GetX)
        .method("SetX", &Body::Track::SetX)
        .method("ReverseTrack", &Body::Track::ReverseTrack);

    mod.add_type<Body>("Body")
        .method("density", &Body::density)
        .method("ye", &Body::ye);

    // ========== Vacuum ==========

    mod.add_type<Vacuum::Track>("VacuumTrack", jlcxx::julia_base_type<Body::Track>())
        .constructor<double>()
        .constructor<double, double>()
        .constructor<double, double, double>();

    // Factory for VacuumTrack
    mod.method("make_VacuumTrack", [](double xend) {
        return std::make_shared<Vacuum::Track>(xend);
    });
    mod.method("make_VacuumTrack_2", [](double xini, double xend) {
        return std::make_shared<Vacuum::Track>(xini, xend);
    });

    mod.add_type<Vacuum>("Vacuum", jlcxx::julia_base_type<Body>())
        .constructor<>();

    // Factory function returning shared_ptr for use with Set_Body
    mod.method("make_Vacuum", []() { return std::make_shared<Vacuum>(); });

    // ========== ConstantDensity ==========

    mod.add_type<ConstantDensity::Track>("ConstantDensityTrack", jlcxx::julia_base_type<Body::Track>())
        .constructor<double>()
        .constructor<double, double>();

    mod.add_type<ConstantDensity>("ConstantDensity", jlcxx::julia_base_type<Body>())
        .constructor<double, double>();

    mod.method("make_ConstantDensity", [](double density, double ye) {
        return std::make_shared<ConstantDensity>(density, ye);
    });

    // Factory for ConstantDensityTrack
    mod.method("make_ConstantDensityTrack", [](double xend) {
        return std::make_shared<ConstantDensity::Track>(xend);
    });
    mod.method("make_ConstantDensityTrack_2", [](double xini, double xend) {
        return std::make_shared<ConstantDensity::Track>(xini, xend);
    });

    // ========== VariableDensity ==========

    mod.add_type<VariableDensity::Track>("VariableDensityTrack", jlcxx::julia_base_type<Body::Track>())
        .constructor<double>()
        .constructor<double, double>();

    mod.add_type<VariableDensity>("VariableDensity", jlcxx::julia_base_type<Body>())
        .constructor<std::vector<double>, std::vector<double>, std::vector<double>>();

    mod.method("make_VariableDensity", [](const std::vector<double>& x,
                                          const std::vector<double>& density,
                                          const std::vector<double>& ye) {
        return std::make_shared<VariableDensity>(x, density, ye);
    });

    // Factory for VariableDensityTrack
    mod.method("make_VariableDensityTrack", [](double xend) {
        return std::make_shared<VariableDensity::Track>(xend);
    });
    mod.method("make_VariableDensityTrack_2", [](double xini, double xend) {
        return std::make_shared<VariableDensity::Track>(xini, xend);
    });

    // ========== Earth ==========

    mod.add_type<Earth::Track>("EarthTrack", jlcxx::julia_base_type<Body::Track>())
        .constructor<double>()  // baseline only
        .constructor<double, double, double>();  // xini, xend, baseline

    mod.add_type<Earth>("Earth", jlcxx::julia_base_type<Body>())
        .constructor<std::string>()
        .constructor<std::vector<double>, std::vector<double>, std::vector<double>>()
        .method("GetRadius", &Earth::GetRadius);

    mod.method("make_Earth", [](const std::string& model) {
        return std::make_shared<Earth>(model);
    });
    mod.method("make_Earth_default", []() {
        return std::make_shared<Earth>();
    });

    // Factory for EarthTrack
    mod.method("make_EarthTrack", [](double baseline) {
        return std::make_shared<Earth::Track>(baseline);
    });
    mod.method("make_EarthTrack_3", [](double xini, double xend, double baseline) {
        return std::make_shared<Earth::Track>(xini, xend, baseline);
    });

    // ========== EarthAtm ==========

    mod.add_type<EarthAtm::Track>("EarthAtmTrack", jlcxx::julia_base_type<Body::Track>())
        .constructor<double, double, double>()
        .constructor<double, double, double, double>();

    mod.add_type<EarthAtm>("EarthAtm", jlcxx::julia_base_type<Body>())
        .constructor<std::string>()
        .constructor<std::vector<double>, std::vector<double>, std::vector<double>>()
        .method("GetRadius", &EarthAtm::GetRadius)
        .method("GetAtmosphereHeight", &EarthAtm::GetAtmosphereHeight)
        .method("SetAtmosphereHeight", &EarthAtm::SetAtmosphereHeight)
        .method("MakeTrack", &EarthAtm::MakeTrack)
        .method("MakeTrackWithCosine", &EarthAtm::MakeTrackWithCosine);

    mod.method("make_EarthAtm", [](const std::string& model) {
        return std::make_shared<EarthAtm>(model);
    });
    mod.method("make_EarthAtm_default", []() {
        return std::make_shared<EarthAtm>();
    });

    // Factory for EarthAtmTrack
    mod.method("make_EarthAtmTrack", [](double phi_nu, double costh_nu, double production_height) {
        return std::make_shared<EarthAtm::Track>(phi_nu, costh_nu, production_height);
    });
    mod.method("make_EarthAtmTrack_4", [](double phi_nu, double costh_nu, double production_height, double atm_height) {
        return std::make_shared<EarthAtm::Track>(phi_nu, costh_nu, production_height, atm_height);
    });

    // ========== Sun ==========

    mod.add_type<Sun::Track>("SunTrack", jlcxx::julia_base_type<Body::Track>())
        .constructor<double>()
        .constructor<double, double>();

    mod.add_type<Sun>("Sun", jlcxx::julia_base_type<Body>())
        .constructor<std::string, bool>();

    mod.method("make_Sun", [](const std::string& model, bool use_composition) {
        return std::make_shared<Sun>(model, use_composition);
    });
    mod.method("make_Sun_default", []() {
        return std::make_shared<Sun>();
    });

    // Factory for SunTrack
    mod.method("make_SunTrack", [](double xend) {
        return std::make_shared<Sun::Track>(xend);
    });
    mod.method("make_SunTrack_2", [](double xini, double xend) {
        return std::make_shared<Sun::Track>(xini, xend);
    });

    // ========== SunASnu ==========

    mod.add_type<SunASnu::Track>("SunASnuTrack", jlcxx::julia_base_type<Body::Track>())
        .constructor<double>()
        .constructor<double, double>();

    mod.add_type<SunASnu>("SunASnu", jlcxx::julia_base_type<Body>())
        .constructor<std::string, bool>();

    mod.method("make_SunASnu", [](const std::string& model, bool use_composition) {
        return std::make_shared<SunASnu>(model, use_composition);
    });
    mod.method("make_SunASnu_default", []() {
        return std::make_shared<SunASnu>();
    });

    // Factory for SunASnuTrack
    mod.method("make_SunASnuTrack", [](double xend) {
        return std::make_shared<SunASnu::Track>(xend);
    });
    mod.method("make_SunASnuTrack_2", [](double xini, double xend) {
        return std::make_shared<SunASnu::Track>(xini, xend);
    });

    // ========== Cross Sections ==========

    mod.add_type<NeutrinoCrossSections>("NeutrinoCrossSections");

    mod.add_type<NullCrossSections>("NullCrossSections", jlcxx::julia_base_type<NeutrinoCrossSections>())
        .constructor<>();

    mod.add_type<GlashowResonanceCrossSection>("GlashowResonanceCrossSection",
        jlcxx::julia_base_type<NeutrinoCrossSections>())
        .constructor<>();

    mod.add_type<NeutrinoDISCrossSectionsFromTables>("NeutrinoDISCrossSectionsFromTables",
        jlcxx::julia_base_type<NeutrinoCrossSections>())
        .constructor<std::string>()
        .method("GetEmin", &NeutrinoDISCrossSectionsFromTables::GetEmin)
        .method("GetEmax", &NeutrinoDISCrossSectionsFromTables::GetEmax);

    mod.add_type<CrossSectionLibrary>("CrossSectionLibrary")
        .constructor<>()
        .method("numberOfTargets", &CrossSectionLibrary::numberOfTargets);

    // Cross section loading functions
    mod.method("loadDefaultCrossSections", &loadDefaultCrossSections);
    mod.method("loadWCG24CrossSections", &loadWCG24CrossSections);
    mod.method("loadWCG24NuclearCrossSections", &loadWCG24NuclearCrossSections);

    // ========== TauDecaySpectra ==========

    mod.add_type<TauDecaySpectra>("TauDecaySpectra")
        .method("GetTauToHadronBranchingRatio", &TauDecaySpectra::GetTauToHadronBranchingRatio)
        .method("GetTauToLeptonBranchingRatio", &TauDecaySpectra::GetTauToLeptonBranchingRatio);

    // ========== nuSQUIDS Main Class ==========

    mod.add_type<nuSQUIDS>("nuSQUIDS")
        // Constructors
        .constructor<>()
        .constructor<unsigned int, NeutrinoType>()
        .constructor<std::string>()
        .constructor<std::string, std::string>()
        // Core methods
        .method("EvolveState", &nuSQUIDS::EvolveState)
        .method("GetNumNeu", &nuSQUIDS::GetNumNeu)
        .method("GetNumRho", &nuSQUIDS::GetNumRho)
        // Set methods - use lambdas to handle shared_ptr conversion for body types
        .method("Set_Track", &nuSQUIDS::Set_Track)
        .method("Set_E", &nuSQUIDS::Set_E)
        .method("Set_rel_error", &nuSQUIDS::Set_rel_error)
        .method("Set_abs_error", &nuSQUIDS::Set_abs_error)
        .method("Set_h_min", &nuSQUIDS::Set_h_min)
        .method("Set_h_max", &nuSQUIDS::Set_h_max)
        .method("Set_h", &nuSQUIDS::Set_h)
        // Get methods
        .method("GetTrack", &nuSQUIDS::GetTrack)
        .method("GetBody", &nuSQUIDS::GetBody)
        // Mixing parameters
        .method("Set_MixingAngle", &nuSQUIDS::Set_MixingAngle)
        .method("Set_SquareMassDifference", &nuSQUIDS::Set_SquareMassDifference)
        .method("Set_CPPhase", &nuSQUIDS::Set_CPPhase)
        .method("Get_MixingAngle", &nuSQUIDS::Get_MixingAngle)
        .method("Get_SquareMassDifference", &nuSQUIDS::Get_SquareMassDifference)
        .method("Get_CPPhase", &nuSQUIDS::Get_CPPhase)
        // HDF5 I/O
        .method("WriteStateHDF5", [](nuSQUIDS& nus, const std::string& filename) {
            nus.WriteStateHDF5(filename);
        })
        .method("ReadStateHDF5", [](nuSQUIDS& nus, const std::string& filename) {
            nus.ReadStateHDF5(filename);
        })
        // Interaction settings
        .method("Set_TauRegeneration", &nuSQUIDS::Set_TauRegeneration)
        .method("Set_GlashowResonance", &nuSQUIDS::Set_GlashowResonance)
        .method("Set_PositivityConstrain", &nuSQUIDS::Set_PositivityConstrain)
        .method("Set_PositivityConstrainStep", &nuSQUIDS::Set_PositivityConstrainStep)
        // Progress bar
        .method("Set_ProgressBar", &nuSQUIDS::Set_ProgressBar);

    // Set_Body wrappers for each body type (handle shared_ptr conversion)
    mod.method("Set_Body_Vacuum", [](nuSQUIDS& nus, std::shared_ptr<Vacuum> body) {
        nus.Set_Body(std::static_pointer_cast<Body>(body));
    });
    mod.method("Set_Body_ConstantDensity", [](nuSQUIDS& nus, std::shared_ptr<ConstantDensity> body) {
        nus.Set_Body(std::static_pointer_cast<Body>(body));
    });
    mod.method("Set_Body_VariableDensity", [](nuSQUIDS& nus, std::shared_ptr<VariableDensity> body) {
        nus.Set_Body(std::static_pointer_cast<Body>(body));
    });
    mod.method("Set_Body_Earth", [](nuSQUIDS& nus, std::shared_ptr<Earth> body) {
        nus.Set_Body(std::static_pointer_cast<Body>(body));
    });
    mod.method("Set_Body_EarthAtm", [](nuSQUIDS& nus, std::shared_ptr<EarthAtm> body) {
        nus.Set_Body(std::static_pointer_cast<Body>(body));
    });
    mod.method("Set_Body_Sun", [](nuSQUIDS& nus, std::shared_ptr<Sun> body) {
        nus.Set_Body(std::static_pointer_cast<Body>(body));
    });
    mod.method("Set_Body_SunASnu", [](nuSQUIDS& nus, std::shared_ptr<SunASnu> body) {
        nus.Set_Body(std::static_pointer_cast<Body>(body));
    });

    // Set_Track wrappers for each track type (handle shared_ptr conversion)
    mod.method("Set_Track_Vacuum", [](nuSQUIDS& nus, std::shared_ptr<Vacuum::Track> track) {
        nus.Set_Track(std::static_pointer_cast<Body::Track>(track));
    });
    mod.method("Set_Track_ConstantDensity", [](nuSQUIDS& nus, std::shared_ptr<ConstantDensity::Track> track) {
        nus.Set_Track(std::static_pointer_cast<Body::Track>(track));
    });
    mod.method("Set_Track_VariableDensity", [](nuSQUIDS& nus, std::shared_ptr<VariableDensity::Track> track) {
        nus.Set_Track(std::static_pointer_cast<Body::Track>(track));
    });
    mod.method("Set_Track_Earth", [](nuSQUIDS& nus, std::shared_ptr<Earth::Track> track) {
        nus.Set_Track(std::static_pointer_cast<Body::Track>(track));
    });
    mod.method("Set_Track_EarthAtm", [](nuSQUIDS& nus, std::shared_ptr<EarthAtm::Track> track) {
        nus.Set_Track(std::static_pointer_cast<Body::Track>(track));
    });
    mod.method("Set_Track_Sun", [](nuSQUIDS& nus, std::shared_ptr<Sun::Track> track) {
        nus.Set_Track(std::static_pointer_cast<Body::Track>(track));
    });
    mod.method("Set_Track_SunASnu", [](nuSQUIDS& nus, std::shared_ptr<SunASnu::Track> track) {
        nus.Set_Track(std::static_pointer_cast<Body::Track>(track));
    });

    // Evaluation methods (need separate registration due to overloads)
    mod.method("EvalFlavor_single", [](const nuSQUIDS& nus, unsigned int flv) {
        return nus.EvalFlavor(flv);
    });
    mod.method("EvalFlavor_energy", [](const nuSQUIDS& nus, unsigned int flv, double E) {
        return nus.EvalFlavor(flv, E);
    });
    mod.method("EvalFlavor_energy_rho", [](const nuSQUIDS& nus, unsigned int flv, double E, unsigned int rho) {
        return nus.EvalFlavor(flv, E, rho);
    });

    mod.method("EvalMass_single", [](const nuSQUIDS& nus, unsigned int mass_state) {
        return nus.EvalMass(mass_state);
    });
    // EvalMass with energy requires rho argument (no 2-arg version exists)
    mod.method("EvalMass_energy_rho", [](const nuSQUIDS& nus, unsigned int mass_state, double E, unsigned int rho) {
        return nus.EvalMass(mass_state, E, rho);
    });

    // Set initial state methods (need wrappers for marray conversion)
    mod.method("Set_initial_state_1D", [](nuSQUIDS& nus, const std::vector<double>& state, Basis basis) {
        nus.Set_initial_state(vec_to_marray1(state), basis);
    });

    mod.method("Set_initial_state_2D", [](nuSQUIDS& nus, const std::vector<double>& state,
                                          size_t dim1, size_t dim2, Basis basis) {
        nus.Set_initial_state(vec_to_marray2(state, dim1, dim2), basis);
    });

    mod.method("Set_initial_state_3D", [](nuSQUIDS& nus, const std::vector<double>& state,
                                          size_t dim1, size_t dim2, size_t dim3, Basis basis) {
        nus.Set_initial_state(vec_to_marray3(state, dim1, dim2, dim3), basis);
    });

    // Constructor with energy array (uses helper for marray conversion)
    mod.method("nuSQUIDS_multi", [](const std::vector<double>& E_vec, unsigned int numneu,
                                    NeutrinoType NT, bool iinteraction) {
        return std::make_shared<nuSQUIDS>(vec_to_marray1(E_vec), numneu, NT, iinteraction);
    });

    mod.method("nuSQUIDS_multi_noint", [](const std::vector<double>& E_vec, unsigned int numneu,
                                          NeutrinoType NT) {
        return std::make_shared<nuSQUIDS>(vec_to_marray1(E_vec), numneu, NT);
    });

    mod.method("nuSQUIDS_multi_xsec", [](const std::vector<double>& E_vec, unsigned int numneu,
                                         NeutrinoType NT, bool iinteraction,
                                         std::shared_ptr<CrossSectionLibrary> ncs) {
        return std::make_shared<nuSQUIDS>(vec_to_marray1(E_vec), numneu, NT, iinteraction, ncs);
    });

    // Get energy range as vector
    mod.method("GetERange", [](const nuSQUIDS& nus) -> std::vector<double> {
        auto arr = nus.GetERange();
        std::vector<double> result(arr.size());
        std::copy(arr.get_data(), arr.get_data() + arr.size(), result.begin());
        return result;
    });

    // GSL step function setter
    mod.method("Set_GSL_step", [](nuSQUIDS& nus, GSL_STEP_FUNCTIONS step) {
        wrap_Set_GSL_STEP(&nus, step);
    });

    // ========== nuSQUIDSAtm Class ==========

    mod.add_type<nuSQUIDSAtm<nuSQUIDS>>("nuSQUIDSAtm")
        .method("EvolveState", &nuSQUIDSAtm<nuSQUIDS>::EvolveState)
        .method("GetNumNeu", &nuSQUIDSAtm<nuSQUIDS>::GetNumNeu)
        // Use lambdas to resolve overloaded methods
        .method("Set_rel_error", [](nuSQUIDSAtm<nuSQUIDS>& atm, double eps) { atm.Set_rel_error(eps); })
        .method("Set_abs_error", [](nuSQUIDSAtm<nuSQUIDS>& atm, double eps) { atm.Set_abs_error(eps); })
        .method("Set_TauRegeneration", &nuSQUIDSAtm<nuSQUIDS>::Set_TauRegeneration)
        .method("Set_GlashowResonance", &nuSQUIDSAtm<nuSQUIDS>::Set_GlashowResonance)
        .method("Set_ProgressBar", &nuSQUIDSAtm<nuSQUIDS>::Set_ProgressBar)
        .method("GetNumCos", &nuSQUIDSAtm<nuSQUIDS>::GetNumCos)
        .method("GetNumE", &nuSQUIDSAtm<nuSQUIDS>::GetNumE);

    // nuSQUIDSAtm constructors
    mod.method("nuSQUIDSAtm_new", [](const std::vector<double>& costh_vec,
                                     const std::vector<double>& E_vec,
                                     unsigned int numneu, NeutrinoType NT,
                                     bool iinteraction) {
        return std::make_shared<nuSQUIDSAtm<nuSQUIDS>>(
            vec_to_marray1(costh_vec), vec_to_marray1(E_vec), numneu, NT, iinteraction);
    });

    mod.method("nuSQUIDSAtm_new_xsec", [](const std::vector<double>& costh_vec,
                                          const std::vector<double>& E_vec,
                                          unsigned int numneu, NeutrinoType NT,
                                          bool iinteraction,
                                          std::shared_ptr<CrossSectionLibrary> ncs) {
        return std::make_shared<nuSQUIDSAtm<nuSQUIDS>>(
            vec_to_marray1(costh_vec), vec_to_marray1(E_vec), numneu, NT, iinteraction, ncs);
    });

    // nuSQUIDSAtm evaluation
    mod.method("EvalFlavor_atm", [](const nuSQUIDSAtm<nuSQUIDS>& atm, unsigned int flv,
                                    double costh, double E, unsigned int rho) {
        return atm.EvalFlavor(flv, costh, E, rho);
    });

    // nuSQUIDSAtm set initial state
    mod.method("Set_initial_state_atm", [](nuSQUIDSAtm<nuSQUIDS>& atm,
                                           const std::vector<double>& state,
                                           size_t dim1, size_t dim2, size_t dim3,
                                           Basis basis) {
        atm.Set_initial_state(vec_to_marray3(state, dim1, dim2, dim3), basis);
    });

    mod.method("Set_initial_state_atm_4D", [](nuSQUIDSAtm<nuSQUIDS>& atm,
                                              const std::vector<double>& state,
                                              size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                                              Basis basis) {
        marray<double, 4> arr;
        arr.resize(std::array<size_t,4>{dim1, dim2, dim3, dim4});
        std::copy(state.begin(), state.end(), arr.get_data());
        atm.Set_initial_state(arr, basis);
    });

    // ========== Utility Functions ==========

    mod.method("getResourcePath", &getResourcePath);
    mod.method("linspace", &linspace_vec);
    mod.method("logspace", &logspace_vec);

    // PDG code helpers
    mod.method("isNuclearPDGCode", &isNuclearPDGCode);
    mod.method("getAtomicNumber", &getAtomicNumber);
    mod.method("getMassNumber", &getMassNumber);
}
