"""
    nuSQuIDS

Julia bindings for the nuSQuIDS neutrino evolution library.

nuSQuIDS solves the neutrino evolution equations in various environments
(vacuum, constant density, Earth, Sun, etc.) and optionally includes
non-coherent interactions (CC, NC, tau regeneration, Glashow resonance).

# Quick Start

```julia
using nuSQuIDS

# Create physical constants
units = Const()

# Single energy mode - oscillation probability
nus = nuSQUIDS(3, neutrino)
Set_Body(nus, Vacuum())
Set_Track(nus, VacuumTrack(100.0 * km(units)))
Set_E(nus, 1.0 * GeV(units))
Set_initial_state(nus, [0.0, 1.0, 0.0], flavor)  # pure nu_mu
EvolveState(nus)
prob = EvalFlavor(nus, 0)  # P(nu_mu -> nu_e)
println("P(nu_mu -> nu_e) = \$prob")
```

# Multiple Energy Mode

```julia
using nuSQuIDS

units = Const()
E_nodes = logspace(1e9, 1e12, 100)  # 1 GeV to 1 TeV
nus = nuSQUIDS_multi(E_nodes, 3, both, true)  # with interactions
# ... set body, track, initial flux, then evolve
```

See the examples/ directory for more comprehensive usage examples.
"""
module nuSQuIDS

using CxxWrap
using CxxWrap.StdLib: StdVector
using Libdl

# Path to the compiled C++ library
function get_lib_path()
    # Try different possible locations
    candidates = [
        joinpath(@__DIR__, "..", "..", "lib", "nuSQuIDSjl"),
        joinpath(@__DIR__, "..", "..", "..", "lib", "nuSQuIDSjl"),
        "/usr/local/lib/nuSQuIDSjl"
    ]

    for path in candidates
        for ext in [".dylib", ".so", ""]
            full_path = path * ext
            if isfile(full_path)
                return dirname(full_path), basename(path)
            end
        end
    end

    error("Could not find nuSQuIDSjl library. Please build the Julia bindings first.")
end

const _lib_dir, _lib_name = get_lib_path()

# Load the C++ wrapper module
@wrapmodule(() -> joinpath(_lib_dir, _lib_name), :define_julia_module, Libdl.RTLD_GLOBAL)

function __init__()
    @initcxx
end

# ============== Exports ==============

# Enums
export Basis, mass, flavor, interaction
export NeutrinoType, neutrino, antineutrino, both
export GSLStepFunction, GSL_STEP_RK2, GSL_STEP_RK4, GSL_STEP_RKF45
export GSL_STEP_RKCK, GSL_STEP_RK8PD, GSL_STEP_MSADAMS
export NeutrinoFlavor, electron_flavor, muon_flavor, tau_flavor, sterile_flavor
export Current, CC, NC, GR
export PDGCode

# Physical constants
export Const, GF, Na, PeV, TeV, GeV, MeV, keV, eV
export kg, gr, meter, cm, km, fermi, angstrom, AU, parsec
export picobarn, femtobarn, sec, hour, day, year

# SU_vector
export SU_vector, Dim, Size, SetAllComponents, GetComponents
export Transpose, Imag, Real, Projector, PosProjector, NegProjector
export SU_vector_Identity, SU_vector_Generator

# Bodies
export Body, BodyTrack
export Vacuum, VacuumTrack, make_Vacuum, make_VacuumTrack, make_VacuumTrack_2
export ConstantDensity, ConstantDensityTrack, make_ConstantDensity, make_ConstantDensityTrack, make_ConstantDensityTrack_2
export VariableDensity, VariableDensityTrack, make_VariableDensity, make_VariableDensityTrack, make_VariableDensityTrack_2
export Earth, EarthTrack, make_Earth, make_Earth_default, make_EarthTrack, make_EarthTrack_3
export EarthAtm, EarthAtmTrack, MakeTrack, MakeTrackWithCosine, make_EarthAtm, make_EarthAtm_default
export make_EarthAtmTrack, make_EarthAtmTrack_4
export Sun, SunTrack, make_Sun, make_Sun_default, make_SunTrack, make_SunTrack_2
export SunASnu, SunASnuTrack, make_SunASnu, make_SunASnu_default, make_SunASnuTrack, make_SunASnuTrack_2
export GetRadius, GetAtmosphereHeight, SetAtmosphereHeight
export Set_Track_Vacuum, Set_Track_ConstantDensity, Set_Track_VariableDensity
export Set_Track_Earth, Set_Track_EarthAtm, Set_Track_Sun, Set_Track_SunASnu

# Track methods
export GetInitialX, GetFinalX, GetX, SetX, ReverseTrack

# Cross sections
export NeutrinoCrossSections, NullCrossSections
export GlashowResonanceCrossSection, NeutrinoDISCrossSectionsFromTables
export CrossSectionLibrary, numberOfTargets
export loadDefaultCrossSections, loadWCG24CrossSections, loadWCG24NuclearCrossSections
export GetEmin, GetEmax

# TauDecaySpectra
export TauDecaySpectra, GetTauToHadronBranchingRatio, GetTauToLeptonBranchingRatio

# nuSQUIDS main class
export nuSQUIDS, nuSQUIDS_multi, nuSQUIDS_multi_noint, nuSQUIDS_multi_xsec
export EvolveState, GetNumNeu, GetNumRho
export Set_Body, Set_Track, Set_E
export Set_Body_Vacuum, Set_Body_ConstantDensity, Set_Body_VariableDensity
export Set_Body_Earth, Set_Body_EarthAtm, Set_Body_Sun, Set_Body_SunASnu
export Set_rel_error, Set_abs_error, Set_h_min, Set_h_max, Set_h
export GetTrack, GetBody
export Set_MixingAngle, Set_SquareMassDifference, Set_CPPhase
export Get_MixingAngle, Get_SquareMassDifference, Get_CPPhase
export WriteStateHDF5, ReadStateHDF5
export Set_TauRegeneration, Set_GlashowResonance
export Set_PositivityConstrain, Set_PositivityConstrainStep
export Set_ProgressBar, Set_GSL_step
export GetERange

# Evaluation methods
export EvalFlavor, EvalMass
export EvalFlavor_single, EvalFlavor_energy, EvalFlavor_energy_rho
export EvalMass_single, EvalMass_energy, EvalMass_energy_rho

# Set initial state
export Set_initial_state, Set_initial_state_1D, Set_initial_state_2D, Set_initial_state_3D

# nuSQUIDSAtm
export nuSQUIDSAtm, nuSQUIDSAtm_new, nuSQUIDSAtm_new_xsec
export GetNumCos, GetNumE
export EvalFlavor_atm
export Set_initial_state_atm, Set_initial_state_atm_4D

# Utilities
export getResourcePath, linspace, logspace
export isNuclearPDGCode, getAtomicNumber, getMassNumber
export density, ye

# ============== Helper Functions ==============

"""
    to_std_vector(arr::AbstractVector{Float64}) -> StdVector{Float64}

Convert a Julia array to a C++ std::vector for passing to C++ functions.
"""
function to_std_vector(arr::AbstractVector{Float64})
    sv = StdVector{Float64}()
    for x in arr
        push!(sv, x)
    end
    return sv
end

# ============== Convenience Functions ==============

"""
    EvalFlavor(nus::nuSQUIDS, flavor::Integer)
    EvalFlavor(nus::nuSQUIDS, flavor::Integer, E::Real)
    EvalFlavor(nus::nuSQUIDS, flavor::Integer, E::Real, rho::Integer)

Evaluate the flavor composition of the neutrino state.

In single energy mode, use the first form.
In multiple energy mode, specify the energy E (in eV).
For neutrino+antineutrino mode, specify rho (0=neutrino, 1=antineutrino).
"""
function EvalFlavor(nus::nuSQUIDS, flavor::Integer)
    return EvalFlavor_single(nus, UInt32(flavor))
end

function EvalFlavor(nus::nuSQUIDS, flavor::Integer, E::Real)
    return EvalFlavor_energy(nus, UInt32(flavor), Float64(E))
end

function EvalFlavor(nus::nuSQUIDS, flavor::Integer, E::Real, rho::Integer)
    return EvalFlavor_energy_rho(nus, UInt32(flavor), Float64(E), UInt32(rho))
end

"""
    EvalMass(nus::nuSQUIDS, mass_state::Integer)
    EvalMass(nus::nuSQUIDS, mass_state::Integer, E::Real)
    EvalMass(nus::nuSQUIDS, mass_state::Integer, E::Real, rho::Integer)

Evaluate the mass eigenstate composition.
"""
function EvalMass(nus::nuSQUIDS, mass_state::Integer)
    return EvalMass_single(nus, UInt32(mass_state))
end

function EvalMass(nus::nuSQUIDS, mass_state::Integer, E::Real)
    return EvalMass_energy(nus, UInt32(mass_state), Float64(E))
end

function EvalMass(nus::nuSQUIDS, mass_state::Integer, E::Real, rho::Integer)
    return EvalMass_energy_rho(nus, UInt32(mass_state), Float64(E), UInt32(rho))
end

"""
    Set_initial_state(nus::nuSQUIDS, state::Vector{Float64}, basis::Basis)

Set the initial neutrino state (single energy mode).
"""
function Set_initial_state(nus::nuSQUIDS, state::Vector{Float64}, basis)
    Set_initial_state_1D(nus, to_std_vector(state), basis)
end

"""
    Set_initial_state(nus::nuSQUIDS, state::Matrix{Float64}, basis::Basis)

Set the initial neutrino state (multiple energy mode, neutrino or antineutrino only).
State shape: (n_energies, n_flavors)
"""
function Set_initial_state(nus::nuSQUIDS, state::Matrix{Float64}, basis)
    dim1, dim2 = size(state)
    # Julia is column-major, C++ is row-major - flatten appropriately
    Set_initial_state_2D(nus, to_std_vector(vec(transpose(state))), dim1, dim2, basis)
end

"""
    Set_initial_state(nus::nuSQUIDS, state::Array{Float64,3}, basis::Basis)

Set the initial neutrino state (multiple energy mode with neutrino+antineutrino).
State shape: (n_energies, 2, n_flavors)
"""
function Set_initial_state(nus::nuSQUIDS, state::Array{Float64,3}, basis)
    dim1, dim2, dim3 = size(state)
    # Flatten for C++ (row-major order)
    flat = Float64[]
    for i in 1:dim1
        for j in 1:dim2
            for k in 1:dim3
                push!(flat, state[i, j, k])
            end
        end
    end
    Set_initial_state_3D(nus, to_std_vector(flat), dim1, dim2, dim3, basis)
end

"""
    EvalFlavor(atm::nuSQUIDSAtm, flavor::Integer, costh::Real, E::Real, rho::Integer)

Evaluate the flavor composition for atmospheric neutrinos.
"""
function EvalFlavor(atm, flavor::Integer, costh::Real, E::Real, rho::Integer)
    return EvalFlavor_atm(atm, UInt32(flavor), Float64(costh), Float64(E), UInt32(rho))
end

"""
    Set_initial_state(atm::nuSQUIDSAtm, state::Array{Float64,3}, basis::Basis)

Set the initial state for atmospheric neutrinos.
State shape: (n_costh, n_energies, n_flavors)
"""
function Set_initial_state(atm, state::Array{Float64,3}, basis)
    dim1, dim2, dim3 = size(state)
    flat = Float64[]
    for i in 1:dim1
        for j in 1:dim2
            for k in 1:dim3
                push!(flat, state[i, j, k])
            end
        end
    end
    Set_initial_state_atm(atm, to_std_vector(flat), dim1, dim2, dim3, basis)
end

"""
    Set_initial_state(atm::nuSQUIDSAtm, state::Array{Float64,4}, basis::Basis)

Set the initial state for atmospheric neutrinos (neutrino+antineutrino mode).
State shape: (n_costh, n_energies, 2, n_flavors)
"""
function Set_initial_state(atm, state::Array{Float64,4}, basis)
    dim1, dim2, dim3, dim4 = size(state)
    flat = Float64[]
    for i in 1:dim1
        for j in 1:dim2
            for k in 1:dim3
                for l in 1:dim4
                    push!(flat, state[i, j, k, l])
                end
            end
        end
    end
    Set_initial_state_atm_4D(atm, to_std_vector(flat), dim1, dim2, dim3, dim4, basis)
end

end # module
