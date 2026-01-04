# Astrophysical Data Files

This directory contains data files for astrophysical body models used by nuSQuIDS.

## Earth Model Files

### EARTH_MODEL_PREM.dat (Isoscalar Mode)

The Preliminary Reference Earth Model (PREM) in isoscalar format.

**Format:** 3 columns, space-separated
| Column | Description | Units |
|--------|-------------|-------|
| 1 | Relative radius (r/R) | dimensionless [0, 1] |
| 2 | Density | g/cm^3 |
| 3 | Electron fraction (Ye) | dimensionless |

The relative radius runs from 0 (center of Earth) to 1 (surface).

### EARTH_MODEL_PREM_wIso.dat (Nuclear Composition Mode)

The PREM model with detailed nuclear composition for each element.

**Format:** 13 columns, space-separated
| Column | Description | Units |
|--------|-------------|-------|
| 1 | Relative radius (r/R) | dimensionless [0, 1] |
| 2 | Density | g/cm^3 |
| 3 | Electron fraction (Ye) | dimensionless |
| 4 | Hydrogen (H) fraction | dimensionless |
| 5 | Oxygen (O) fraction | dimensionless |
| 6 | Sodium (Na) fraction | dimensionless |
| 7 | Magnesium (Mg) fraction | dimensionless |
| 8 | Aluminum (Al) fraction | dimensionless |
| 9 | Silicon (Si) fraction | dimensionless |
| 10 | Sulfur (S) fraction | dimensionless |
| 11 | Calcium (Ca) fraction | dimensionless |
| 12 | Iron (Fe) fraction | dimensionless |
| 13 | Nickel (Ni) fraction | dimensionless |

The element fractions represent number fractions (not mass fractions) and should sum to approximately 1 for each row.

## Sun Model Files

### bs05_agsop.dat, bs05op.dat

Standard Solar Model (SSM) files from Bahcall, Serenelli & Basu (2005).

**Reference:** Bahcall, Serenelli, Basu, "New Solar Opacities, Abundances, Helioseismology, and Neutrino Fluxes", Astrophys.J. 621:L85-L88 (2005), [arXiv:astro-ph/0412440](https://arxiv.org/abs/astro-ph/0412440)

- **bs05op.dat**: Standard solar model with older (OP) opacities
- **bs05_agsop.dat**: Standard solar model with AGS 2005 metallicity and OP opacities

**Format:** 12 columns, space-separated

| Column | Description | Units |
|--------|-------------|-------|
| 1 | Mass fraction | M/M_sun [0, 1] |
| 2 | Radius | R/R_sun [0, 1] |
| 3 | Temperature | K |
| 4 | Density | g/cm^3 |
| 5 | Pressure | dyn/cm^2 |
| 6 | Luminosity fraction | L/L_sun [0, 1] |
| 7 | Hydrogen (H) | mass fraction |
| 8 | Helium-4 (He4) | mass fraction |
| 9 | Helium-3 (He3) | mass fraction |
| 10 | Carbon-12 (C12) | mass fraction |
| 11 | Nitrogen-14 (N14) | mass fraction |
| 12 | Oxygen-16 (O16) | mass fraction |

The mass fractions of elements (columns 7-12) represent the dominant species in the solar interior that are relevant for nuclear reactions and opacity calculations.

### bs05op-org.dat

Original unmodified version of the bs05op.dat file.

### nele_bs05op.dat

Electron number density profile for the Standard Solar Model.

**Format:** 2 columns, space-separated

| Column | Description | Units |
|--------|-------------|-------|
| 1 | Radius | R/R_sun [0, 1] |
| 2 | Electron number density | cm^-3 |

This file is used by the `Sun` and `SunASnu` body classes for calculating the matter potential experienced by neutrinos propagating through the solar interior.

## Usage

When creating custom body models, users can provide either:

1. **Isoscalar mode**: 3 columns (r/R, density, Ye)
   - Used with standard neutrino cross sections on isoscalar nucleon targets

2. **Nuclear composition mode**: 13 columns (r/R, density, Ye, H, O, Na, Mg, Al, Si, S, Ca, Fe, Ni)
   - Used with nuclear-specific cross sections for detailed composition tracking

The `Earth` and `EarthAtm` classes automatically detect the file format based on the number of columns.
