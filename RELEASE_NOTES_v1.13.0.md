# nuSQuIDS v1.13.0 Release Notes

**Release Date:** January 3, 2026

This is a major feature release that introduces nuclear composition support for precision neutrino propagation calculations, along with significant improvements to documentation, installation, and Python bindings.

## Acknowledgments

We thank **Philip Weigel** and **Alex Wen** for implementing and reviewing the initial nuclear cross section and composition infrastructure that this release builds upon.

## Major New Features

### Nuclear Composition Support

This release introduces comprehensive support for nuclear composition in neutrino propagation calculations, enabling precision studies of matter effects with per-nucleus cross sections.

**Body Composition:**
- `Earth` and `EarthAtm` classes now support detailed nuclear composition (H, O, Na, Mg, Al, Si, S, Ca, Fe, Ni)
- New PREM Earth model file with isotopic composition: `EARTH_MODEL_PREM_wIso.dat`
- Automatic detection of file format (3-column isoscalar vs 13-column composition mode)
- `composition()` method on Body classes returns element number fractions at any position

**Sun Composition:**
- `Sun` and `SunASnu` classes now support composition via `use_composition_information` parameter
- Parses Standard Solar Model files for H, He4, He3, C12, N14, O16 mass fractions
- Automatic conversion from mass fractions to number fractions
- Based on Bahcall, Serenelli & Basu (2005) SSM data

**New PDG Codes:**
- `helium3` (1000020030)
- `helium4` (1000020040)
- `nitrogen` (1000070140)

### WCG24 Nuclear Cross Sections

Complete set of per-nucleus neutrino DIS cross sections from Weigel, Conrad & Garcia (2024) [arXiv:2408.05866]:

| Target | File |
|--------|------|
| Isoscalar | `wcg24_base_isoscalar.h5` |
| Proton | `wcg24_base_proton.h5` |
| Neutron | `wcg24_base_neutron.h5` |
| Oxygen-16 | `wcg24_oxygen.h5` |
| Carbon-12 | `wcg24_carbon.h5` |
| Sodium-23 | `wcg24_sodium.h5` |
| Magnesium-24 | `wcg24_magnesium.h5` |
| Aluminum-27 | `wcg24_aluminum.h5` |
| Silicon-28 | `wcg24_silicon.h5` |
| Sulfur-32 | `wcg24_sulfur.h5` |
| Calcium-40 | `wcg24_calcium.h5` |
| Iron-56 | `wcg24_iron.h5` |
| Nickel-58 | `wcg24_nickel.h5` |
| Lead-208 | `wcg24_lead.h5` |

New helper functions:
- `loadWCG24CrossSections()` - Load WCG24 proton/neutron cross sections
- `loadWCG24NuclearCrossSections()` - Load all nuclear cross sections for Earth composition

## Bug Fixes

- **Glashow resonance electron density:** Fixed incorrect electron density calculation when using nuclear cross sections. The Glashow resonance now correctly uses electron number density derived from composition.
- **GSL 2.8 compatibility:** Fixed compilation issues with GSL 2.8.
- **Linux build:** Removed macOS-only linker flags that caused build failures on Linux.
- **Python linking:** Fixed `python-config --ldflags` usage on Linux for proper linking.
- **Binder compatibility:** Fixed postBuild script to handle empty CONDA_PREFIX.

## Documentation Improvements

- **New README files:**
  - `examples/README.md` - Comprehensive guide to all 15 example programs
  - `benchmarks/README.md` - Benchmark suite documentation
  - `data/xsections/README.md` - Cross section data file documentation
  - `data/astro/README.md` - Earth and Sun model file documentation

- **Sun model documentation:** Added detailed format specification for SSM files (bs05_agsop.dat, bs05op.dat) with reference to Bahcall et al. (2005).

- **Main README updates:** Added citation section with BibTeX entries for nuSQuIDS and SQuIDS papers.

## Installation Improvements

- **pip installation:** Full support for `pip install .` with automatic SQuIDS installation
- **Package name:** Changed pip package name from 'nusquids' to 'nuSQuIDS' for consistency
- **pybind11 bindings:** Now the default for Python bindings (Boost.Python still supported)
- **Binder support:** Updated to use pybind11, pinned Python version for stability

## Python Bindings

- **Comprehensive docstrings:** All major classes and methods now have Python docstrings
- **Python examples:** Added Python versions for all compatible C++ examples
- **New class alias:** `nuSQuIDS` class alias for `nuSQUIDS` (lowercase compatibility)
- **Composition support:** Full Python bindings for composition methods and nuclear cross sections

## New Tests

- `nuclear_xs` - Nuclear cross section integration test
- `xs_regression` - Cross section regression test with validated baselines from v1.12.2 and nuclearxs-dev

## New Examples

- `Composition/` - Body composition for nuclear cross sections
- `EarthComposition/` - Isoscalar vs nuclear composition comparison

## Benchmark Suite

New benchmark infrastructure (`make benchmark` / `make benchmark-quick`) testing:
- Single energy mode
- Multiple energy mode (with/without interactions)
- Atmospheric mode (with/without interactions, Glashow, tau regeneration)
- Nuclear composition mode with PREM

## API Changes

### New Methods
- `Body::composition(const GenericTrack&)` - Returns element number fractions
- `CrossSectionLibrary::numberOfTargets()` - Get target count
- `CrossSectionLibrary::targets()` - Get list of target PDG codes

### New Constructors
- `Sun(std::string sunmodel, bool use_composition_information = false)`
- `SunASnu(std::string sunmodel, bool use_composition_information = false)`
- `ConstantDensity(double density, double ye, std::map<PDGCode, double> composition)`
- `VariableDensity(std::vector<double> x, std::vector<double> density, std::vector<double> ye, std::map<PDGCode, std::vector<double>> composition)`

### Backward Compatibility
All changes are backward compatible. Existing code will continue to work without modification.

## Building

```bash
# C++ library
./configure
make
make test
make install

# Python bindings
pip install .
```

## Citation

If you use nuSQuIDS in your research, please cite:

```bibtex
@article{Arguelles:2021twb,
    author = {Argüelles, Carlos A. and Salvado, Jordi and Weaver, Christopher N.},
    title = "{nuSQuIDS: A toolbox for neutrino propagation}",
    doi = "10.1016/j.cpc.2022.108346",
    journal = "Comput. Phys. Commun.",
    volume = "277",
    pages = "108346",
    year = "2022"
}
```

For nuclear cross sections, please also cite:

```bibtex
@article{Weigel:2024gzh,
    author = "Weigel, Philip L. R. and Conrad, Janet M. and Garcia-Soto, Alfonso",
    title = "{Cross sections and inelasticity distributions of high-energy neutrino deep inelastic scattering}",
    doi = "10.1103/PhysRevD.111.043044",
    journal = "Phys. Rev. D",
    volume = "111",
    pages = "043044",
    year = "2025"
}
```
