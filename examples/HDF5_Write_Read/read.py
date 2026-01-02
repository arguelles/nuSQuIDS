#!/usr/bin/env python3
"""
HDF5 Read Example
=================

This example demonstrates how to read nuSQUIDS state from HDF5 files
and use the information to extract fluxes.

This is the Python equivalent of read.cpp in this directory.
Run write.py first to generate the HDF5 files.
"""

import numpy as np
import nuSQuIDS as nsq


def main():
    units = nsq.Const()

    # Read states from HDF5 files
    print("Reading initial state from initial_state.hdf5...")
    inus = nsq.nuSQUIDS("./initial_state.hdf5")

    print("Reading final state from final_state.hdf5...")
    fnus = nsq.nuSQUIDS("./final_state.hdf5")

    # Get number of neutrinos
    numneu = inus.GetNumNeu()
    print(f"Number of neutrino flavors: {numneu}")

    # Write flux information to text file
    print("\nWriting fluxes to fluxes_flavor.txt...")

    Nen = 1000
    lEmin = 0.0  # log10(E/GeV)
    lEmax = 4.0

    with open("fluxes_flavor.txt", "w") as f:
        # Write header
        header = "# log10(E/GeV) E[eV]"
        for fl in range(numneu):
            header += f" flux_{fl} ratio_{fl}"
        f.write(header + "\n")

        for lE in np.linspace(lEmin, lEmax, Nen):
            E = 10**lE * units.GeV
            line = f"{lE:.6f} {E:.6e}"

            for fl in range(numneu):
                # Get final flux
                flux_final = fnus.EvalFlavor(fl, E, 0)
                # Get initial muon flux (index 1) for ratio calculation
                flux_initial = inus.EvalFlavor(1, E, 0)

                if flux_initial > 0:
                    ratio = flux_final / flux_initial
                else:
                    ratio = 0.0

                line += f" {flux_final:.6e} {ratio:.6e}"

            f.write(line + "\n")

    print("Done! Results written to fluxes_flavor.txt")

    # Optionally plot results
    try:
        import matplotlib.pyplot as plt

        data = np.loadtxt("fluxes_flavor.txt")
        lE = data[:, 0]

        plt.figure(figsize=(10, 6))

        flavor_names = [r'$\nu_e$', r'$\nu_\mu$', r'$\nu_\tau$']
        colors = ['blue', 'red', 'green']

        for i, (name, color) in enumerate(zip(flavor_names, colors)):
            ratio_col = 3 + 2 * i  # ratio columns are at indices 3, 5, 7
            plt.plot(lE, data[:, ratio_col], label=name, color=color)

        plt.xlabel(r'$\log_{10}(E/\mathrm{GeV})$')
        plt.ylabel(r'Flux ratio (final/initial $\nu_\mu$)')
        plt.title('Neutrino oscillation through Earth (upgoing)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig('hdf5_read_plot.png', dpi=150, bbox_inches='tight')
        plt.show()
        print("Plot saved to hdf5_read_plot.png")

    except ImportError:
        print("Matplotlib not available. Install with: pip install matplotlib")


if __name__ == "__main__":
    main()
