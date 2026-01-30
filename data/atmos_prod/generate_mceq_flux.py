#!/usr/bin/env python3
"""
Generate atmospheric neutrino differential flux production profiles using MCEq.

This script produces the data file needed by EmittingEarthAtm in nuSQuIDS.
Output format: coszen height energy nu_e nu_e_bar nu_mu nu_mu_bar

Requirements:
    pip install MCEq crflux

Usage:
    python generate_mceq_flux.py
"""

from MCEq.core import config, MCEqRun
import crflux.models as crf
import numpy as np
import sys

# Configuration
OUTPUT_FILE = "PROD_MODEL_MCEQ.dat"
MONTH = 'January'
LOCATION = 'SouthPole'

# Grid parameters
H_PTS = 100          # Number of height points
COSZEN_PTS = 50      # Number of cosine zenith points
H_MIN_CM = 1e2       # 1 meter in cm
H_MAX_CM = 6e6       # 60 km in cm

def main():
    print(f"Generating atmospheric flux data...")
    print(f"  Height points: {H_PTS}")
    print(f"  Coszen points: {COSZEN_PTS}")
    print(f"  Location: {LOCATION}, {MONTH}")

    # Create grids
    coszen_grid = np.linspace(-1, 1, COSZEN_PTS)
    zen_grid = np.arccos(coszen_grid)
    zen_deg_grid = np.degrees(zen_grid)

    # Height grid from 60km to 1m (in cm), logarithmically spaced
    h_grid = np.logspace(np.log10(H_MAX_CM), np.log10(H_MIN_CM), H_PTS)

    # Initialize MCEq
    print("Initializing MCEq...")
    mceq = MCEqRun(
        interaction_model='SIBYLL23C',
        primary_model=(crf.HillasGaisser2012, 'H3a'),
        theta_deg=zen_deg_grid[0],
        density_model=('MSIS00_IC', (LOCATION, MONTH)),
    )

    # Open output file
    with open(OUTPUT_FILE, 'w') as f:
        # Write header comment
        f.write("# Atmospheric neutrino differential flux production from MCEq\n")
        f.write("# Format: coszen height[cm] energy[GeV] nu_e nu_e_bar nu_mu nu_mu_bar\n")
        f.write(f"# Generated with: SIBYLL23C, HillasGaisser2012 H3a, {LOCATION} {MONTH}\n")
        f.write(f"# Grid: {COSZEN_PTS} coszen x {H_PTS} heights x {len(mceq.e_grid)} energies\n")

        total_iterations = COSZEN_PTS * H_PTS
        current = 0

        for izen, coszen in enumerate(coszen_grid):
            # Set zenith angle
            mceq.set_theta_deg(zen_deg_grid[izen])

            # Convert heights to slant depths
            X_grid = mceq.density_model.h2X(h_grid)

            # Solve cascade equations
            mceq.solve(int_grid=X_grid)

            for idx in range(H_PTS):
                # Get fluxes for all flavors (neutrinos and anti-neutrinos)
                nue_flux = mceq.get_solution('nue', grid_idx=idx)
                nuebar_flux = mceq.get_solution('antinue', grid_idx=idx)
                numu_flux = mceq.get_solution('numu', grid_idx=idx)
                numubar_flux = mceq.get_solution('antinumu', grid_idx=idx)

                for ien, energy in enumerate(mceq.e_grid):
                    # Write: coszen height energy nu_e nu_e_bar nu_mu nu_mu_bar
                    f.write(f"{coszen:.6e} {h_grid[idx]:.6e} {energy:.6e} "
                           f"{nue_flux[ien]:.6e} {nuebar_flux[ien]:.6e} "
                           f"{numu_flux[ien]:.6e} {numubar_flux[ien]:.6e}\n")

                current += 1
                if current % 100 == 0:
                    print(f"  Progress: {current}/{total_iterations} ({100*current/total_iterations:.1f}%)")

    print(f"Done! Output written to {OUTPUT_FILE}")
    print(f"File contains {COSZEN_PTS * H_PTS * len(mceq.e_grid)} data points")

if __name__ == "__main__":
    main()
