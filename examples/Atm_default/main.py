#!/usr/bin/env python3
"""
Atmospheric Neutrino Example
============================

This example demonstrates how to use the nuSQUIDSAtm class to solve
the atmospheric neutrino propagation problem, which is both energy
and zenith-angle dependent.

The nuSQUIDSAtm class solves a set of nuSQUIDS energy-dependent objects
for different zenith angles to handle the full 2D problem.

This is the Python equivalent of main.cpp in this directory.
"""

import numpy as np
import nuSQuIDS as nsq


def flux_function(enu, cz):
    """Initial flux function. Returns 1.0 for all energies and angles."""
    return 1.0


def main():
    units = nsq.Const()

    # Get user input for number of neutrino flavors
    print("Enter number of neutrino flavors to consider:")
    print("  (3) Three Active Neutrinos")
    print("  (4) 3+1 Three Active and One Sterile Neutrino")

    numneu = int(input("Enter 3 or 4: "))
    if numneu not in [3, 4]:
        raise ValueError("Only 3 or 4 are valid options")

    # Ask about interactions
    use_int = input("Consider Earth absorption and interactions? (yes/no): ")
    interactions = use_int.lower() in ['yes', 'y']

    # Energy and zenith ranges
    Emin = 1.0e1 * units.GeV
    Emax = 1.0e6 * units.GeV
    czmin = -1.0
    czmax = 0.0

    # Create nuSQUIDSAtm object
    print("\nConstructing nuSQUIDS-Atm object...")
    cz_nodes = nsq.linspace(czmin, czmax, 40)
    E_nodes = nsq.logspace(Emin, Emax, 100)

    nus_atm = nsq.nuSQUIDSAtm(cz_nodes, E_nodes, numneu, nsq.NeutrinoType.both, interactions)
    print("Done constructing nuSQUIDS-Atm object")

    # Set mixing parameters
    print("\nSetting mixing angles...")
    nus_atm.Set_MixingAngle(0, 1, 0.563942)  # theta_12
    nus_atm.Set_MixingAngle(0, 2, 0.154085)  # theta_13
    nus_atm.Set_MixingAngle(1, 2, 0.785398)  # theta_23

    nus_atm.Set_SquareMassDifference(1, 7.65e-05)  # dm^2_21
    nus_atm.Set_SquareMassDifference(2, 0.00247)   # dm^2_31

    nus_atm.Set_CPPhase(0, 2, 0.0)

    if numneu > 3:
        # Sterile neutrino parameters
        nus_atm.Set_SquareMassDifference(3, -1.0)  # dm^2_41
        nus_atm.Set_MixingAngle(1, 3, 0.160875)    # theta_24

    print("Done setting mixing angles")

    # Setup integration precision
    nus_atm.Set_rel_error(1.0e-6)
    nus_atm.Set_abs_error(1.0e-6)
    nus_atm.Set_GSL_step(nsq.GSL_STEP_RK4)

    # Get energy and zenith ranges
    e_range = nus_atm.GetERange()
    cz_range = nus_atm.GetCosthRange()

    # Construct the initial state
    print("\nSetting initial state...")

    # Initial state array: [cos_zenith, energy, neutrino/antineutrino, flavor]
    num_cz = nus_atm.GetNumCos()
    num_E = nus_atm.GetNumE()

    # Create 4D array for initial state
    inistate = np.zeros((num_cz, num_E, 2, numneu))

    for ci in range(num_cz):
        for ei in range(num_E):
            for rho in range(2):  # 0=neutrino, 1=antineutrino
                for flv in range(numneu):
                    # Set muon neutrino flux only (flv=1)
                    if flv == 1:
                        inistate[ci, ei, rho, flv] = flux_function(e_range[ei], cz_range[ci])
                    else:
                        inistate[ci, ei, rho, flv] = 0.0

    # Set the initial state
    nus_atm.Set_initial_state(inistate, nsq.Basis.flavor)
    print("Done setting initial state")

    # Enable progress bar and oscillations
    nus_atm.Set_ProgressBar(True)
    nus_atm.Set_IncludeOscillations(True)

    # Evolve the state
    print("\nEvolving state...")
    nus_atm.EvolveState()
    print("Done evolving")

    # Write results to file
    print("\nWriting results to fluxes_flavor.txt...")
    Nen = 700
    Ncz = 100
    lEmin = np.log10(Emin / units.GeV)
    lEmax = np.log10(Emax / units.GeV)

    with open("fluxes_flavor.txt", "w") as f:
        # Write header
        header = "# log10(E/GeV) cos(zenith) E[GeV]"
        for fl in range(numneu):
            header += f" nu_{fl}"
        for fl in range(numneu):
            header += f" nubar_{fl}"
        f.write(header + "\n")

        for cz in np.linspace(czmin, czmax, Ncz):
            for lE in np.linspace(lEmin, lEmax, Nen):
                E = 10**lE * units.GeV
                line = f"{lE:.6f} {cz:.6f} {10**lE:.6e}"

                # Neutrino flavors
                for fl in range(numneu):
                    flux = nus_atm.EvalFlavor(fl, cz, E, 0)
                    line += f" {flux:.6e}"

                # Antineutrino flavors
                for fl in range(numneu):
                    flux = nus_atm.EvalFlavor(fl, cz, E, 1)
                    line += f" {flux:.6e}"

                f.write(line + "\n")
            f.write("\n")  # Blank line between zenith angles (for gnuplot)

    print("Done! Results written to fluxes_flavor.txt")

    # Optionally plot results
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        # Read and reshape data for plotting
        # Plot muon neutrino survival probability at a fixed zenith angle

        print("\nCreating plots...")

        # Sample some zenith angles and energies
        cz_values = [-0.9, -0.5, -0.1]
        lE_array = np.linspace(lEmin, lEmax, 200)

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Plot neutrinos
        ax = axes[0]
        for cz in cz_values:
            mu_flux = []
            for lE in lE_array:
                E = 10**lE * units.GeV
                mu_flux.append(nus_atm.EvalFlavor(1, cz, E, 0))  # mu neutrino

            ax.plot(lE_array, mu_flux, label=f'cos(z) = {cz:.1f}')

        ax.set_xlabel(r'$\log_{10}(E/\mathrm{GeV})$')
        ax.set_ylabel(r'$\nu_\mu$ flux ratio')
        ax.set_title(r'Muon neutrino survival')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot antineutrinos
        ax = axes[1]
        for cz in cz_values:
            mu_flux = []
            for lE in lE_array:
                E = 10**lE * units.GeV
                mu_flux.append(nus_atm.EvalFlavor(1, cz, E, 1))  # mu antineutrino

            ax.plot(lE_array, mu_flux, label=f'cos(z) = {cz:.1f}')

        ax.set_xlabel(r'$\log_{10}(E/\mathrm{GeV})$')
        ax.set_ylabel(r'$\bar{\nu}_\mu$ flux ratio')
        ax.set_title(r'Muon antineutrino survival')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('atmospheric_oscillations.png', dpi=150, bbox_inches='tight')
        plt.show()
        print("Plot saved to atmospheric_oscillations.png")

    except ImportError:
        print("Matplotlib not available. Install with: pip install matplotlib")


if __name__ == "__main__":
    main()
