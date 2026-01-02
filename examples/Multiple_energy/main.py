#!/usr/bin/env python3
"""
Multiple Energy Mode Example
============================

This example demonstrates how to use nuSQUIDS to propagate neutrinos
considering oscillations for the multiple energy case. The neutrinos
are propagated through the Earth with a power-law initial spectrum.

This is the Python equivalent of main.cpp in this directory.
"""

import numpy as np
import nuSQuIDS as nsq

def main():
    units = nsq.Const()

    # Number of neutrinos: 3 for standard, 4 for 3+1 sterile
    print("Select number of neutrino flavors:")
    print("  (3) Three Active Neutrinos")
    print("  (4) 3+1 Three Active and One Sterile Neutrino")

    numneu = int(input("Enter 3 or 4: "))
    if numneu not in [3, 4]:
        raise ValueError("Only 3 or 4 are valid options")

    # Create energy array (logarithmically spaced from 1 GeV to 10 TeV)
    E_min = 1.0 * units.GeV
    E_max = 1.0e4 * units.GeV
    E_nodes = nsq.logspace(E_min, E_max, 200)

    # Create nuSQUIDS object
    # Arguments: energy array, number of neutrinos, neutrino type, interactions
    nus = nsq.nuSQUIDS(E_nodes, numneu, nsq.NeutrinoType.neutrino, False)

    # Set up the trajectory through the Earth
    # Zenith angle: cos(theta) = -1 means upgoing (through the Earth)
    phi = np.arccos(-1.0)

    # Create Earth atmospheric model and track
    earth_atm = nsq.EarthAtm()
    track_atm = earth_atm.MakeTrack(phi)

    nus.Set_Body(earth_atm)
    nus.Set_Track(track_atm)

    # Set mixing angles (in radians)
    nus.Set_MixingAngle(0, 1, 0.563942)  # theta_12
    nus.Set_MixingAngle(0, 2, 0.154085)  # theta_13
    nus.Set_MixingAngle(1, 2, 0.785398)  # theta_23

    # Set mass squared differences (in eV^2)
    nus.Set_SquareMassDifference(1, 7.65e-05)  # dm^2_21
    nus.Set_SquareMassDifference(2, 0.00247)   # dm^2_31

    if numneu == 4:
        # Sterile neutrino parameters
        nus.Set_SquareMassDifference(3, 0.1)  # dm^2_41
        nus.Set_MixingAngle(1, 3, 0.1)        # theta_24

    # Set maximum step size for integration
    nus.Set_h_max(500.0 * units.km)

    # Set GSL step function (RK4)
    nus.Set_GSL_step(nsq.GSL_STEP_RK4)

    # Set numerical precision
    nus.Set_rel_error(1.0e-5)
    nus.Set_abs_error(1.0e-5)

    # Enable progress bar
    nus.Set_ProgressBar(True)

    # Construct the initial state
    E_range = nus.GetERange()

    # Create initial state array: power-law spectrum for muon neutrinos
    N0 = 1.0e18
    inistate = np.zeros((len(E_range), numneu))
    for i, E in enumerate(E_range):
        inistate[i, 1] = N0 * E**(-2)  # muon neutrino flux (index 1)

    # Set initial state
    nus.Set_initial_state(inistate, nsq.Basis.flavor)

    # Propagate the neutrinos through the Earth
    print("\nEvolving neutrino state...")
    nus.EvolveState()

    print("\nWriting outputs...")

    # Save results to file
    Nen = 1000
    lEmin = 0.0  # log10(E/GeV)
    lEmax = 4.0

    with open("fluxes_flavor.txt", "w") as f:
        header = "# log10(E/GeV)  E[eV]  " + "  ".join(
            [f"flux_{i}  ratio_{i}" for i in range(numneu)]
        )
        f.write(header + "\n")

        for lE in np.linspace(lEmin, lEmax, Nen):
            E = 10**lE * units.GeV
            line = f"{lE:.6f} {E:.6e}"

            for fl in range(numneu):
                flux = nus.EvalFlavor(fl, E)
                initial_flux = N0 * E**(-2)
                ratio = flux / initial_flux if initial_flux > 0 else 0
                line += f" {flux:.6e} {ratio:.6e}"

            f.write(line + "\n")

    print("\nDone! Results written to fluxes_flavor.txt")

    # Optionally create a simple plot
    try:
        import matplotlib.pyplot as plt

        # Read back the data for plotting
        data = np.loadtxt("fluxes_flavor.txt")
        lE = data[:, 0]

        plt.figure(figsize=(10, 6))

        flavor_names = ['e', 'mu', 'tau'] + (['s'] if numneu == 4 else [])
        colors = ['blue', 'red', 'green', 'orange']

        for i, (name, color) in enumerate(zip(flavor_names, colors)):
            ratio_col = 3 + 2*i  # ratio columns are at indices 3, 5, 7, ...
            plt.plot(lE, data[:, ratio_col], label=rf'$\nu_{{{name}}}$', color=color)

        plt.xlabel(r'$\log_{10}(E/\mathrm{GeV})$')
        plt.ylabel(r'Flux ratio (final/initial)')
        plt.title('Neutrino oscillation through Earth (upgoing)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig('oscillation_plot.png', dpi=150, bbox_inches='tight')
        plt.show()
        print("Plot saved to oscillation_plot.png")

    except ImportError:
        print("Matplotlib not available. Install with: pip install matplotlib")


if __name__ == "__main__":
    main()
