#!/usr/bin/env python3
"""
HDF5 Write Example
==================

This example demonstrates how to use nuSQUIDS to save the full state
of the system to HDF5 files. The state is saved both before and after
propagation through the Earth.

This is the Python equivalent of write.cpp in this directory.
"""

import numpy as np
import nuSQuIDS as nsq


def main():
    units = nsq.Const()

    # Number of neutrinos
    numneu = 3

    # Create nuSQUIDS object for multiple energy mode
    # 200 energy bins from 1 GeV to 10 TeV
    E_nodes = nsq.logspace(1.0 * units.GeV, 1.0e4 * units.GeV, 200)
    nus = nsq.nuSQUIDS(E_nodes, numneu, nsq.NeutrinoType.neutrino, False)

    # Set up trajectory through the Earth
    # zenith angle: cos(theta) = -1 means straight up through Earth
    phi = np.arccos(-1.0)

    earth_atm = nsq.EarthAtm()
    track_atm = earth_atm.MakeTrack(phi)

    nus.Set_Body(earth_atm)
    nus.Set_Track(track_atm)

    # Set mixing parameters
    nus.Set_MixingAngle(0, 1, 0.563942)  # theta_12
    nus.Set_MixingAngle(0, 2, 0.154085)  # theta_13
    nus.Set_MixingAngle(1, 2, 0.785398)  # theta_23
    nus.Set_SquareMassDifference(1, 7.65e-05)  # dm^2_21
    nus.Set_SquareMassDifference(2, 0.00247)   # dm^2_31

    # Set integration parameters
    nus.Set_h_max(500.0 * units.km)
    nus.Set_GSL_step(nsq.GSL_STEP_RK4)
    nus.Set_rel_error(1.0e-5)
    nus.Set_abs_error(1.0e-5)

    # Enable progress bar
    nus.Set_ProgressBar(True)

    # Construct the initial state
    E_range = nus.GetERange()
    N0 = 1.0e18

    # Power-law spectrum for muon neutrinos
    inistate = np.zeros((len(E_range), numneu))
    for i, E in enumerate(E_range):
        inistate[i, 1] = N0 * E**(-2)  # muon neutrino (index 1)

    # Set initial state
    nus.Set_initial_state(inistate, nsq.Basis.flavor)

    # Save initial state to HDF5
    print("Writing initial state to initial_state.hdf5...")
    nus.WriteStateHDF5("./initial_state.hdf5")

    # Propagate through the Earth
    print("\nEvolving state...")
    nus.EvolveState()

    # Save final state to HDF5
    print("\nWriting final state to final_state.hdf5...")
    nus.WriteStateHDF5("./final_state.hdf5")

    print("\nDone! State files written successfully.")


if __name__ == "__main__":
    main()
