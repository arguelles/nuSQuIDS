#!/usr/bin/env python3
"""
Generate atmospheric neutrino differential production profiles using MCEq.

This script produces the data file needed by EmittingEarthAtm in nuSQuIDS.
It computes the derivative of the flux with respect to height to get the
production rate at each point in the atmosphere.

Output format: coszen height[cm] energy[GeV] nu_e nu_e_bar nu_mu nu_mu_bar
Where the flux values are differential production rates (d(flux)/d(height)).

Requirements:
    pip install MCEq crflux numpy

Usage:
    python generate_production_profile.py [--output OUTPUT_FILE]
"""

from MCEq.core import config, MCEqRun
import crflux.models as crf
import numpy as np
import argparse
import sys

# Default configuration
DEFAULT_OUTPUT = "PROD_MODEL_MCEQ.dat"
MONTH = 'January'
LOCATION = 'SouthPole'

# Grid parameters
H_PTS = 100          # Number of height points
COSZEN_PTS = 50      # Number of cosine zenith points
H_MIN_CM = 1e2       # 1 meter in cm
H_MAX_CM = 6e6       # 60 km in cm

def compute_differential(flux_array, height_array):
    """
    Compute the differential production rate from integrated flux.

    Uses central differences for interior points and forward/backward
    differences at boundaries.

    Parameters:
        flux_array: 1D array of flux values at each height
        height_array: 1D array of heights (in cm)

    Returns:
        1D array of differential production rates
    """
    n = len(flux_array)
    diff = np.zeros(n)

    # Use log-spacing aware differentiation
    log_h = np.log(height_array)

    for i in range(n):
        if i == 0:
            # Forward difference at top of atmosphere
            dflux = flux_array[i+1] - flux_array[i]
            dlogh = log_h[i+1] - log_h[i]
            # d(flux)/d(h) = d(flux)/d(log h) * d(log h)/d(h) = d(flux)/d(log h) / h
            diff[i] = dflux / (dlogh * height_array[i]) if dlogh != 0 else 0
        elif i == n - 1:
            # Backward difference at ground
            dflux = flux_array[i] - flux_array[i-1]
            dlogh = log_h[i] - log_h[i-1]
            diff[i] = dflux / (dlogh * height_array[i]) if dlogh != 0 else 0
        else:
            # Central difference for interior points
            dflux = flux_array[i+1] - flux_array[i-1]
            dlogh = log_h[i+1] - log_h[i-1]
            diff[i] = dflux / (dlogh * height_array[i]) if dlogh != 0 else 0

    return diff

def main():
    parser = argparse.ArgumentParser(description='Generate MCEq atmospheric production profile')
    parser.add_argument('--output', '-o', default=DEFAULT_OUTPUT,
                        help=f'Output filename (default: {DEFAULT_OUTPUT})')
    parser.add_argument('--coszen-pts', type=int, default=COSZEN_PTS,
                        help=f'Number of cosine zenith points (default: {COSZEN_PTS})')
    parser.add_argument('--height-pts', type=int, default=H_PTS,
                        help=f'Number of height points (default: {H_PTS})')
    args = parser.parse_args()

    print(f"Generating atmospheric neutrino production profile...")
    print(f"  Height points: {args.height_pts}")
    print(f"  Coszen points: {args.coszen_pts}")
    print(f"  Location: {LOCATION}, {MONTH}")
    print(f"  Output: {args.output}")

    # Create grids
    coszen_grid = np.linspace(-1, 1, args.coszen_pts)
    zen_grid = np.arccos(coszen_grid)
    zen_deg_grid = np.degrees(zen_grid)

    # Height grid from 60km to ground (in cm), logarithmically spaced
    # Going from high altitude to low altitude
    h_grid = np.logspace(np.log10(H_MAX_CM), np.log10(H_MIN_CM), args.height_pts)

    # Initialize MCEq
    print("Initializing MCEq...")
    mceq = MCEqRun(
        interaction_model='SIBYLL23C',
        primary_model=(crf.HillasGaisser2012, 'H3a'),
        theta_deg=zen_deg_grid[0],
        density_model=('MSIS00_IC', (LOCATION, MONTH)),
    )

    e_grid = mceq.e_grid
    n_energies = len(e_grid)

    print(f"  Energy points: {n_energies} ({e_grid[0]:.2e} - {e_grid[-1]:.2e} GeV)")

    # Storage for all data
    all_data = []

    total_zenith = args.coszen_pts
    for izen, coszen in enumerate(coszen_grid):
        print(f"  Processing zenith {izen+1}/{total_zenith} (cos(zen) = {coszen:.3f})...")

        # Set zenith angle
        mceq.set_theta_deg(zen_deg_grid[izen])

        # Convert heights to slant depths
        X_grid = mceq.density_model.h2X(h_grid)

        # Solve cascade equations
        mceq.solve(int_grid=X_grid)

        # Get fluxes at all heights for each flavor
        # Shape: (n_heights, n_energies)
        nue_flux = np.array([mceq.get_solution('nue', grid_idx=idx) for idx in range(args.height_pts)])
        nuebar_flux = np.array([mceq.get_solution('antinue', grid_idx=idx) for idx in range(args.height_pts)])
        numu_flux = np.array([mceq.get_solution('numu', grid_idx=idx) for idx in range(args.height_pts)])
        numubar_flux = np.array([mceq.get_solution('antinumu', grid_idx=idx) for idx in range(args.height_pts)])
        nutau_flux = np.array([mceq.get_solution('nutau', grid_idx=idx) for idx in range(args.height_pts)])
        nutaubar_flux = np.array([mceq.get_solution('antinutau', grid_idx=idx) for idx in range(args.height_pts)])

        # Compute differential production for each energy
        for ie in range(n_energies):
            # Get flux column for this energy
            nue_col = nue_flux[:, ie]
            nuebar_col = nuebar_flux[:, ie]
            numu_col = numu_flux[:, ie]
            numubar_col = numubar_flux[:, ie]
            nutau_col = nutau_flux[:, ie]
            nutaubar_col = nutaubar_flux[:, ie]

            # Compute differentials
            dnue = compute_differential(nue_col, h_grid)
            dnuebar = compute_differential(nuebar_col, h_grid)
            dnumu = compute_differential(numu_col, h_grid)
            dnumubar = compute_differential(numubar_col, h_grid)
            dnutau = compute_differential(nutau_col, h_grid)
            dnutaubar = compute_differential(nutaubar_col, h_grid)

            # Store data for each height
            for ih in range(args.height_pts):
                all_data.append([
                    coszen,
                    h_grid[ih],
                    e_grid[ie],
                    dnue[ih],
                    dnuebar[ih],
                    dnumu[ih],
                    dnumubar[ih],
                    dnutau[ih],
                    dnutaubar[ih]
                ])

    # Write output file
    print(f"Writing {len(all_data)} data points to {args.output}...")
    with open(args.output, 'w') as f:
        # Write header
        f.write("# Atmospheric neutrino differential production profile from MCEq\n")
        f.write("# Format: coszen height[cm] energy[GeV] nu_e nu_e_bar nu_mu nu_mu_bar nu_tau nu_tau_bar\n")
        f.write("# Values are differential production rates: d(flux)/d(height)\n")
        f.write(f"# Generated with: SIBYLL23C, HillasGaisser2012 H3a, {LOCATION} {MONTH}\n")
        f.write(f"# Grid: {args.coszen_pts} coszen x {args.height_pts} heights x {n_energies} energies\n")
        f.write("#\n")

        for row in all_data:
            f.write(f"{row[0]:.6e} {row[1]:.6e} {row[2]:.6e} "
                   f"{row[3]:.6e} {row[4]:.6e} {row[5]:.6e} {row[6]:.6e} "
                   f"{row[7]:.6e} {row[8]:.6e}\n")

    print(f"Done! Output written to {args.output}")

if __name__ == "__main__":
    main()
