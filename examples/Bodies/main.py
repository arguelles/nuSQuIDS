#!/usr/bin/env python3
"""
Bodies Example
==============

This example demonstrates how to calculate neutrino oscillation
probabilities for different initial flavor states, using various
predefined bodies and tracks.

Demonstrates:
- Earth (long baseline)
- EarthAtm (atmospheric neutrinos)
- VariableDensity
- Vacuum
- ConstantDensity

This is the Python equivalent of main.cpp in this directory.
Note: The custom exBody example (EarthMod) requires C++ class extension
and is not included here.
"""

import numpy as np
import nuSQuIDS as nsq


def main():
    # Create a nuSQUIDS object for single energy mode
    # Arguments: number of neutrinos, neutrino type
    nus = nsq.nuSQUIDS(3, nsq.NeutrinoType.neutrino)

    # Set mixing parameters
    nus.Set_MixingAngle(0, 1, 0.563942)  # theta_12
    nus.Set_MixingAngle(0, 2, 0.154085)  # theta_13
    nus.Set_MixingAngle(1, 2, 0.785398)  # theta_23
    nus.Set_SquareMassDifference(1, 7.65e-05)  # dm^2_21
    nus.Set_SquareMassDifference(2, 0.00247)   # dm^2_31
    nus.Set_CPPhase(0, 2, 0.0)

    # Define units
    units = nsq.Const()

    # Set initial energy
    nus.Set_E(10.0 * units.GeV)

    # ==============================================================
    # Example 1: Long baseline oscillation through Earth
    # ==============================================================
    print("=" * 50)
    print("Example 1: Earth Long Baseline Neutrino Oscillations")
    print("=" * 50)

    baseline = 500.0 * units.km

    # Create Earth body and track
    earth = nsq.Earth()
    earth_track = nsq.Earth.Track(0.0, baseline, baseline)

    nus.Set_Body(earth)
    nus.Set_Track(earth_track)

    # Set initial state: pure muon neutrino
    ini_state = np.array([0.0, 1.0, 0.0])
    nus.Set_initial_state(ini_state, nsq.Basis.flavor)

    # Set numerical precision
    nus.Set_h_max(200.0 * units.km)
    nus.Set_rel_error(1.0e-12)
    nus.Set_abs_error(1.0e-12)

    print("\nInitial state:")
    for i, flv in enumerate(['nu_e', 'nu_mu', 'nu_tau']):
        print(f"  {flv}: {nus.EvalFlavor(i):.6f}")

    nus.EvolveState()

    print("\nFinal state:")
    for i, flv in enumerate(['nu_e', 'nu_mu', 'nu_tau']):
        print(f"  {flv}: {nus.EvalFlavor(i):.6f}")

    # ==============================================================
    # Example 2: Atmospheric neutrino oscillation
    # ==============================================================
    print("\n" + "=" * 50)
    print("Example 2: Earth Atmospheric Neutrino Oscillations")
    print("=" * 50)

    phi = np.arccos(-1.0)  # Straight up through Earth

    earth_atm = nsq.EarthAtm()
    earth_atm_track = earth_atm.MakeTrack(phi)

    nus.Set_Body(earth_atm)
    nus.Set_Track(earth_atm_track)

    # Set higher precision for atmospheric
    nus.Set_rel_error(1.0e-20)
    nus.Set_abs_error(1.0e-20)

    # Change energy to 100 GeV
    nus.Set_E(100.0 * units.GeV)

    # Reset initial condition
    nus.Set_initial_state(ini_state, nsq.Basis.flavor)

    print("\nInitial state (E = 100 GeV):")
    for i, flv in enumerate(['nu_e', 'nu_mu', 'nu_tau']):
        print(f"  {flv}: {nus.EvalFlavor(i):.6f}")

    nus.EvolveState()

    print("\nFinal state:")
    for i, flv in enumerate(['nu_e', 'nu_mu', 'nu_tau']):
        print(f"  {flv}: {nus.EvalFlavor(i):.6f}")

    # ==============================================================
    # Example 3: Variable density medium
    # ==============================================================
    print("\n" + "=" * 50)
    print("Example 3: Variable Density Neutrino Oscillations")
    print("=" * 50)

    # Define a variable density environment with arbitrary profile
    N = 40
    size = 1000.0 * units.km

    x_arr = np.array([size * (i / N) for i in range(N)])
    density_arr = np.abs(np.cos(np.arange(N, dtype=float)))
    ye_arr = np.abs(np.sin(np.arange(N, dtype=float)))

    vardens = nsq.VariableDensity(x_arr, density_arr, ye_arr)
    track_vardens = nsq.VariableDensity.Track(0.0, 200.0 * units.km)

    nus.Set_Body(vardens)
    nus.Set_Track(track_vardens)

    # Set energy
    nus.Set_E(100.0 * units.GeV)

    # Set initial state to muon neutrino
    ini_state = np.array([0.0, 1.0, 0.0])
    nus.Set_initial_state(ini_state, nsq.Basis.flavor)

    print("\nInitial state:")
    for i, flv in enumerate(['nu_e', 'nu_mu', 'nu_tau']):
        print(f"  {flv}: {nus.EvalFlavor(i):.6f}")

    nus.EvolveState()

    print("\nFinal state:")
    for i, flv in enumerate(['nu_e', 'nu_mu', 'nu_tau']):
        print(f"  {flv}: {nus.EvalFlavor(i):.6f}")

    # ==============================================================
    # Example 4: Vacuum oscillation
    # ==============================================================
    print("\n" + "=" * 50)
    print("Example 4: Vacuum Neutrino Oscillations")
    print("=" * 50)

    vacuum = nsq.Vacuum()
    baseline_vac = 500.0 * units.km
    track_vac = nsq.Vacuum.Track(baseline_vac)

    nus.Set_Body(vacuum)
    nus.Set_Track(track_vac)

    # Set energy to 150 MeV (lower energy to see oscillations)
    nus.Set_E(150.0 * units.MeV)

    # Set initial state to electron neutrino
    ini_state = np.array([1.0, 0.0, 0.0])
    nus.Set_initial_state(ini_state, nsq.Basis.flavor)

    print("\nInitial state (E = 150 MeV):")
    for i, flv in enumerate(['nu_e', 'nu_mu', 'nu_tau']):
        print(f"  {flv}: {nus.EvalFlavor(i):.6f}")

    nus.EvolveState()

    print("\nFinal state:")
    for i, flv in enumerate(['nu_e', 'nu_mu', 'nu_tau']):
        print(f"  {flv}: {nus.EvalFlavor(i):.6f}")

    # ==============================================================
    # Example 5: Constant density medium
    # ==============================================================
    print("\n" + "=" * 50)
    print("Example 5: Constant Density Neutrino Oscillations")
    print("=" * 50)

    # Define constant density environment
    density = 100.0  # g/cm^3
    ye = 0.3         # electron fraction
    constdens = nsq.ConstantDensity(density, ye)
    baseline_cd = 500.0 * units.km
    track_constdens = nsq.ConstantDensity.Track(0.0, baseline_cd)

    nus.Set_Body(constdens)
    nus.Set_Track(track_constdens)

    # Set energy to 10 MeV
    nus.Set_E(10.0 * units.MeV)

    # Set initial state to electron neutrino
    ini_state = np.array([1.0, 0.0, 0.0])
    nus.Set_initial_state(ini_state, nsq.Basis.flavor)

    print("\nInitial state (E = 10 MeV, density = 100 g/cm^3):")
    for i, flv in enumerate(['nu_e', 'nu_mu', 'nu_tau']):
        print(f"  {flv}: {nus.EvalFlavor(i):.6f}")

    nus.EvolveState()

    print("\nFinal state:")
    for i, flv in enumerate(['nu_e', 'nu_mu', 'nu_tau']):
        print(f"  {flv}: {nus.EvalFlavor(i):.6f}")

    print("\n" + "=" * 50)
    print("Done!")
    print("=" * 50)


if __name__ == "__main__":
    main()
