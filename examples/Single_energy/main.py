#!/usr/bin/env python3
"""
Single Energy Mode Example
==========================

This example demonstrates how to calculate neutrino oscillation
probabilities for different initial flavor states, mixing angles
and bodies in the single energy mode.

This is the Python equivalent of main.cpp in this directory.
"""

import nuSQuIDS as nsq

def main():
    # We must first create a nuSQUIDS object. In order to do this
    # we must specify the number of neutrino flavors and if we are
    # going to consider neutrino or antineutrino oscillations.

    # In this example we set N_neutrino = 3 and Type = neutrino.
    nus = nsq.nuSQUIDS(3, nsq.NeutrinoType.neutrino)

    # We will now set the mixing parameters and square mass
    # differences. The angles are given in radians and the
    # mass differences in eV^2.

    # We use the standard parametrization as described in the
    # documentation.

    # Mixing angles
    nus.Set_MixingAngle(0, 1, 0.563942)  # theta_12
    nus.Set_MixingAngle(0, 2, 0.154085)  # theta_13
    nus.Set_MixingAngle(1, 2, 0.785398)  # theta_23

    # Square mass differences
    nus.Set_SquareMassDifference(1, 7.65e-05)  # dm^2_21
    nus.Set_SquareMassDifference(2, 0.00247)   # dm^2_31

    # CP phase
    nus.Set_CPPhase(0, 2, 0.0)

    # Define a Const object for handling units
    units = nsq.Const()

    # Now we set the neutrino energy which we are interested in.
    # Energies are always given in natural units (eV).
    nus.Set_E(10.0 * units.GeV)

    # To calculate atmospheric neutrino oscillation probabilities
    # we need to specify a body (Earth) and a track (trajectory).

    print("*** Earth Atmospheric Neutrino Oscillations ***")

    # Set up the Earth atmospheric model
    # phi is the zenith angle (pi = straight up through Earth)
    import math
    phi = math.acos(-1.0)  # cos(theta) = -1 means upgoing

    earth_atm = nsq.EarthAtm()
    earth_atm_track = earth_atm.MakeTrack(phi)

    nus.Set_Body(earth_atm)
    nus.Set_Track(earth_atm_track)

    # Set the initial state: pure muon neutrino
    # Array format: [nu_e, nu_mu, nu_tau]
    import numpy as np
    ini_state = np.array([0.0, 1.0, 0.0])
    nus.Set_initial_state(ini_state, nsq.Basis.flavor)

    # Setup integration settings
    nus.Set_rel_error(1.0e-20)
    nus.Set_abs_error(1.0e-20)

    # We can change the energy
    nus.Set_E(100.0 * units.GeV)

    # Reset the initial condition
    nus.Set_initial_state(ini_state, nsq.Basis.flavor)

    # Print initial state
    print("Initial state")
    print(f"  E = {100.0} GeV")
    print(f"  nu_e:   {nus.EvalFlavor(0):.6f}")
    print(f"  nu_mu:  {nus.EvalFlavor(1):.6f}")
    print(f"  nu_tau: {nus.EvalFlavor(2):.6f}")

    # Evolve the state through the Earth
    nus.EvolveState()

    # Output the result
    print("\nFinal state (after propagation through Earth)")
    print(f"  E = {100.0} GeV")
    print(f"  nu_e:   {nus.EvalFlavor(0):.6f}")
    print(f"  nu_mu:  {nus.EvalFlavor(1):.6f}")
    print(f"  nu_tau: {nus.EvalFlavor(2):.6f}")

    # Calculate oscillation probabilities
    print("\nOscillation probabilities P(nu_mu -> nu_X):")
    print(f"  P(nu_mu -> nu_e):   {nus.EvalFlavor(0):.6f}")
    print(f"  P(nu_mu -> nu_mu):  {nus.EvalFlavor(1):.6f}")
    print(f"  P(nu_mu -> nu_tau): {nus.EvalFlavor(2):.6f}")


if __name__ == "__main__":
    main()
