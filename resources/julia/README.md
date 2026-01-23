# nuSQuIDS Julia Bindings

Julia bindings for the nuSQuIDS neutrino evolution library using [CxxWrap.jl](https://github.com/JuliaInterop/CxxWrap.jl).

> **⚠️ Status: Work in Progress**
>
> The Julia bindings are currently under development. Building requires updating
> nuSQuIDS headers for C++17 compatibility (CxxWrap.jl requires C++17+, while
> nuSQuIDS was written for C++11). The following deprecated features need updates:
> - `std::result_of` (use `std::invoke_result_t`)
> - `std::iterator` base class (use explicit type aliases)
> - `std::hash<T>::result_type` (use `size_t`)
>
> Contributions welcome!

## Prerequisites

1. **nuSQuIDS** must be installed (headers and library)
2. **Julia** (1.6 or later)
3. **CxxWrap.jl** package

Install CxxWrap.jl in Julia:
```julia
using Pkg
Pkg.add("CxxWrap")
```

## Building

Run the build script:
```bash
cd resources/julia
./build.sh
```

This will:
1. Find the CxxWrap.jl installation
2. Build the C++ wrapper library using CMake
3. Place the library in `resources/julia/lib/`

## Installation

After building, add the package to your Julia environment:

```julia
using Pkg
Pkg.develop(path="/path/to/nuSQuIDS/resources/julia/nuSQuIDS")
```

## Usage

### Basic Example: Vacuum Oscillations

```julia
using nuSQuIDS

# Physical constants
units = Const()

# Create single-energy nuSQUIDS (3 flavors, neutrino)
nus = nuSQUIDS(3, neutrino)

# Set up vacuum propagation
Set_Body(nus, Vacuum())
Set_Track(nus, VacuumTrack(1000.0 * km(units)))  # 1000 km baseline

# Set neutrino energy
Set_E(nus, 1.0 * GeV(units))

# Set mixing parameters
Set_MixingAngle(nus, 0, 1, 0.5836)  # θ₁₂
Set_MixingAngle(nus, 0, 2, 0.1496)  # θ₁₃
Set_MixingAngle(nus, 1, 2, 0.8587)  # θ₂₃
Set_SquareMassDifference(nus, 1, 7.5e-5 * eV(units)^2)  # Δm²₂₁
Set_SquareMassDifference(nus, 2, 2.5e-3 * eV(units)^2)  # Δm²₃₁

# Set initial state: pure νμ
Set_initial_state(nus, [0.0, 1.0, 0.0], flavor)

# Evolve the state
EvolveState(nus)

# Get oscillation probabilities
P_e = EvalFlavor(nus, 0)    # P(νμ → νe)
P_mu = EvalFlavor(nus, 1)   # P(νμ → νμ)
P_tau = EvalFlavor(nus, 2)  # P(νμ → ντ)

println("P(νμ → νe)  = $P_e")
println("P(νμ → νμ)  = $P_mu")
println("P(νμ → ντ)  = $P_tau")
```

### Multiple Energy Mode

```julia
using nuSQuIDS

units = Const()

# Create energy array (1 GeV to 1 TeV, 100 points)
E_nodes = logspace(1e9, 1e12, 100)

# Create nuSQUIDS with interactions
nus = nuSQUIDS_multi(E_nodes, 3, both, true)

# Set Earth body for atmospheric propagation
Set_Body(nus, Earth())

# ... set track, initial flux, mixing parameters ...

# Evolve
EvolveState(nus)

# Evaluate flux at specific energy
flux_numu = EvalFlavor(nus, 1, 100 * GeV(units), 0)  # νμ at 100 GeV
```

### Atmospheric Neutrinos

```julia
using nuSQuIDS

units = Const()

# Zenith angle and energy arrays
costh_nodes = linspace(-1.0, 0.0, 20)  # upgoing
E_nodes = logspace(1e9, 1e15, 100)

# Create atmospheric nuSQUIDS
atm = nuSQUIDSAtm_new(costh_nodes, E_nodes, 3, both, true)

# Set mixing parameters, initial flux...

# Evolve
EvolveState(atm)

# Evaluate at specific zenith and energy
flux = EvalFlavor(atm, 1, -0.5, 1e11, 0)  # νμ at cos(θ)=-0.5, 100 GeV
```

## API Reference

### Enumerations

- `Basis`: `mass`, `flavor`, `interaction`
- `NeutrinoType`: `neutrino`, `antineutrino`, `both`
- `GSLStepFunction`: `GSL_STEP_RK2`, `GSL_STEP_RK4`, `GSL_STEP_RKF45`, etc.
- `Current`: `CC`, `NC`, `GR`

### Bodies

- `Vacuum` - Vacuum propagation
- `ConstantDensity(density, ye)` - Constant matter density
- `VariableDensity(x, density, ye)` - Variable density profile
- `Earth(model_path)` - Earth with PREM model
- `EarthAtm(model_path)` - Earth for atmospheric neutrinos
- `Sun(model_path, use_composition)` - Solar density profile
- `SunASnu(model_path, use_composition)` - Sun for solar atmospheric neutrinos

### Cross Sections

- `loadDefaultCrossSections()` - CSMS cross sections
- `loadWCG24CrossSections()` - WCG24 cross sections (proton/neutron)
- `loadWCG24NuclearCrossSections()` - WCG24 nuclear cross sections

### Utility Functions

- `linspace(min, max, n)` - Linear spacing
- `logspace(min, max, n)` - Logarithmic spacing
- `getResourcePath()` - Path to nuSQuIDS data files

## Troubleshooting

### Library not found

If you get an error about the library not being found, make sure:
1. The build completed successfully
2. The library is in `resources/julia/lib/`
3. nuSQuIDS is installed system-wide

### Symbol errors

If you get symbol not found errors, ensure:
1. nuSQuIDS library version matches the headers used during build
2. Run `sudo make install` to update the system library

## License

This software is licensed under the GNU General Public License v3.0.
