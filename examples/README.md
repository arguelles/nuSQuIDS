# nuSQuIDS Examples

This directory contains example programs demonstrating various features of nuSQuIDS. Most examples are available in both C++ and Python.

## Running Examples

### C++ Examples

After building nuSQuIDS with `make examples`:

```bash
./examples/Single_energy/single_energy
```

### Python Examples

After installing nuSQuIDS (`pip install .`):

```bash
python examples/Single_energy/main.py
```

## Example Overview

| Example | Description | C++ | Python |
|---------|-------------|:---:|:------:|
| [Single_energy](Single_energy/) | Single energy mode oscillations in vacuum | Yes | Yes |
| [Multiple_energy](Multiple_energy/) | Multi-energy propagation with interactions | Yes | Yes |
| [Bodies](Bodies/) | Various body types (Vacuum, Earth, Sun, etc.) | Yes | Yes |
| [Constant_density_layers](Constant_density_layers/) | Multi-layer constant density propagation | Yes | Yes |
| [Atm_default](Atm_default/) | Atmospheric neutrinos using nuSQUIDSAtm | Yes | Yes |
| [Atm_BSM](Atm_BSM/) | BSM physics (NSI) in atmospheric framework | Yes | No |
| [HDF5_Write_Read](HDF5_Write_Read/) | State serialization to HDF5 files | Yes | Yes |
| [Astrophysical_neutrino_flavor_ratio](Astrophysical_neutrino_flavor_ratio/) | Flavor ratio calculation for astrophysical sources | Yes | Yes |
| [NSI](NSI/) | Non-standard neutrino interactions | Yes | No |
| [LV](LV/) | Lorentz violation effects | Yes | No |
| [Xsections](Xsections/) | Custom cross-section implementation | Yes | No |
| [Decoherence](Decoherence/) | Quantum decoherence effects | Yes | No |
| [Composition](Composition/) | Body composition for nuclear cross sections | Yes | Yes |
| [EarthComposition](EarthComposition/) | Isoscalar vs nuclear composition comparison | Yes | Yes |
| [Extended_Source](Extended_Source/) | Extended neutrino source modeling | Yes | No |

## Example Descriptions

### Single_energy

Demonstrates basic neutrino oscillation probability calculations in single energy mode. Shows how to:
- Create a nuSQUIDS object for single energy
- Set mixing parameters and mass-squared differences
- Propagate through vacuum
- Extract oscillation probabilities

### Multiple_energy

Shows multi-energy propagation with a power-law spectrum including:
- Coherent oscillations
- Charged-current (CC) and neutral-current (NC) interactions
- Tau regeneration
- Support for sterile neutrinos (3+1 model)

### Bodies

Demonstrates the various body types available in nuSQuIDS:
- `Vacuum` - Vacuum propagation
- `ConstantDensity` - Uniform matter density
- `VariableDensity` - User-defined density profile
- `Earth` - PREM Earth model (radial propagation)
- `EarthAtm` - Earth for atmospheric neutrinos (chord propagation)
- `Sun` - Standard Solar Model
- `SunASnu` - Sun for solar neutrinos (production to surface)

### Constant_density_layers

Shows propagation through multiple layers of constant density matter, useful for modeling:
- Simplified Earth models
- Detector materials
- Step-function density profiles

### Atm_default

Demonstrates the `nuSQUIDSAtm` class for atmospheric neutrino propagation:
- 2D grid in energy and zenith angle
- Automatic Earth matter effects
- Efficient computation for atmospheric flux calculations

### Atm_BSM

Extends atmospheric neutrino propagation with beyond-Standard-Model physics:
- Non-standard interactions (NSI) in the atmospheric sector
- Custom Hamiltonian modifications

### HDF5_Write_Read

Shows how to save and load nuSQuIDS state to/from HDF5 files:
- Serialization of full propagation state
- Reading saved states for analysis
- Useful for checkpointing long calculations

### Astrophysical_neutrino_flavor_ratio

Calculates flavor ratios for high-energy astrophysical neutrinos:
- Pion decay source (1:2:0)
- Averaged oscillations over cosmological baselines
- Comparison with IceCube measurements

### NSI (C++ only)

Demonstrates non-standard neutrino interactions:
- Matter NSI parameters (epsilon)
- Modified effective potential
- Extending the Hamiltonian with custom operators

### LV (C++ only)

Shows Lorentz violation effects in neutrino propagation:
- CPT-even and CPT-odd operators
- Energy-dependent modifications
- SME (Standard Model Extension) parameters

### Xsections (C++ only)

Demonstrates custom cross-section implementation:
- Subclassing cross-section classes
- Modifying interaction rates
- Testing new physics models

### Decoherence (C++ only)

Models quantum decoherence effects:
- Open quantum system evolution
- Lindblad operators
- Energy-dependent decoherence

### Composition

Shows how to use body composition for nuclear cross sections:
- Setting element abundances (H, O, Fe, etc.)
- Composition-dependent interaction rates
- Relevant for precision calculations

### EarthComposition

Compares isoscalar vs full nuclear composition:
- Default isoscalar approximation
- Full PREM composition with elements
- Quantifying composition effects

### Extended_Source (C++ only)

Models extended neutrino sources:
- Non-point-like emission regions
- Averaging over source extent
- Relevant for nearby sources

## Building Your Own Programs

Use `nusquids-config` to get the correct compiler and linker flags:

```bash
# Compile a C++ program
g++ $(nusquids-config --cflags) myprogram.cpp $(nusquids-config --libs) -o myprogram
```

For Python, simply import nuSQuIDS after installation:

```python
import nuSQuIDS as nsq
```
