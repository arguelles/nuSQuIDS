#!/usr/bin/env python3
"""
Constant Density Layers Example
===============================

This example demonstrates how to propagate neutrinos through multiple
layers with different constant densities. The neutrinos pass through:
1. A vacuum layer (100 km)
2. A constant density layer (50 km, 3.5 g/cm^3, ye=0.5)
3. Another constant density layer (200 km, 10 g/cm^3, ye=0.1)

This is the Python equivalent of main.cpp in this directory.
"""

import numpy as np
import nuSQuIDS as nsq


def main():
    units = nsq.Const()

    # Create nuSQUIDS object for multiple energy mode
    # 60 energy bins from 1 GeV to 10 GeV
    E_nodes = nsq.linspace(1.0 * units.GeV, 10.0 * units.GeV, 60)
    nus = nsq.nuSQUIDS(E_nodes, 3, nsq.NeutrinoType.neutrino, False)

    # Get energy range
    energy_range = nus.GetERange()

    # Set mixing parameters to default
    nus.Set_MixingParametersToDefault()

    # Define the three layers
    # Layer 1: Vacuum, 100 km
    layer_1 = 100.0 * units.km
    vacuum = nsq.Vacuum()
    track_env0 = nsq.Vacuum.Track(layer_1)

    # Layer 2: Constant density, 50 km, 3.5 g/cm^3, ye=0.5
    layer_2 = 50.0 * units.km
    constdens_env1 = nsq.ConstantDensity(3.5, 0.5)
    track_env1 = nsq.ConstantDensity.Track(layer_2)

    # Layer 3: Constant density, 200 km, 10 g/cm^3, ye=0.1
    layer_3 = 200.0 * units.km
    constdens_env2 = nsq.ConstantDensity(10.0, 0.1)
    track_env2 = nsq.ConstantDensity.Track(layer_3)

    # Set the first layer
    nus.Set_Body(vacuum)
    nus.Set_Track(track_env0)

    # Construct the initial state: pure muon neutrino at all energies
    numneu = nus.GetNumNeu()
    numE = nus.GetNumE()
    inistate = np.zeros((numE, numneu))
    for i in range(numE):
        inistate[i, 1] = 1.0  # muon neutrino

    # Set initial state
    nus.Set_initial_state(inistate, nsq.Basis.flavor)

    print("Propagating through Layer 1: Vacuum (100 km)...")
    nus.EvolveState()

    # Switch to second layer and continue evolution
    nus.Set_Body(constdens_env1)
    nus.Set_Track(track_env1)

    print("Propagating through Layer 2: Constant density (50 km, 3.5 g/cm^3)...")
    nus.EvolveState()

    # Switch to third layer and continue evolution
    nus.Set_Body(constdens_env2)
    nus.Set_Track(track_env2)

    print("Propagating through Layer 3: Constant density (200 km, 10 g/cm^3)...")
    nus.EvolveState()

    # Write results to file
    print("\nWriting results to fluxes_flavor.txt...")
    with open("fluxes_flavor.txt", "w") as f:
        f.write("# E[GeV] P(nu_e) P(nu_mu) P(nu_tau)\n")
        for i in range(numE):
            E = energy_range[i]
            line = f"{E / units.GeV:.6f}"
            for k in range(numneu):
                p = nus.EvalFlavorAtNode(k, i)
                line += f" {p:.6e}"
            f.write(line + "\n")

    print("Done!")

    # Optionally plot results
    try:
        import matplotlib.pyplot as plt

        data = np.loadtxt("fluxes_flavor.txt")
        E = data[:, 0]

        plt.figure(figsize=(10, 6))
        flavor_names = [r'$\nu_e$', r'$\nu_\mu$', r'$\nu_\tau$']
        colors = ['blue', 'red', 'green']

        for i, (name, color) in enumerate(zip(flavor_names, colors)):
            plt.plot(E, data[:, i + 1], label=name, color=color)

        plt.xlabel(r'Energy [GeV]')
        plt.ylabel(r'Oscillation probability')
        plt.title('Neutrino oscillation through constant density layers')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig('constant_density_layers.png', dpi=150, bbox_inches='tight')
        plt.show()
        print("Plot saved to constant_density_layers.png")

    except ImportError:
        print("Matplotlib not available. Install with: pip install matplotlib")


if __name__ == "__main__":
    main()
