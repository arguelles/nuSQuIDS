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

### bs05_agsop.dat, bs05op.dat, bs05op-org.dat

Standard Solar Model (SSM) files from Bahcall & Serenelli (2005).

**Format:** Multiple columns including radius, temperature, density, electron density, and nuclear abundances.

### nele_bs05op.dat

Electron number density profile for the Standard Solar Model.

## Usage

When creating custom body models, users can provide either:

1. **Isoscalar mode**: 3 columns (r/R, density, Ye)
   - Used with standard neutrino cross sections on isoscalar nucleon targets

2. **Nuclear composition mode**: 13 columns (r/R, density, Ye, H, O, Na, Mg, Al, Si, S, Ca, Fe, Ni)
   - Used with nuclear-specific cross sections for detailed composition tracking

The `Earth` and `EarthAtm` classes automatically detect the file format based on the number of columns.
