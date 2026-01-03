#!/usr/bin/env python3
"""
nuSQuIDS Body Composition Example

This example demonstrates the use of body composition for nuclear cross sections.

Bodies can provide composition information (elemental fractions) at each point
along a neutrino trajectory. This is useful for:
- Accurate cross section calculations with per-nucleus cross sections
- Studies of neutrino interactions in heterogeneous media
- Precise Earth matter effects with detailed PREM composition

This example:
1. Loads the PREM model with full nuclear composition
2. Demonstrates accessing composition at various Earth depths
3. Shows how to create a ConstantDensity body with custom composition
"""

import nuSQuIDS as nsq
import numpy as np

def get_element_name(code):
    """Get element name from PDGCode."""
    element_names = {
        nsq.hydrogen: "H",
        nsq.oxygen: "O",
        nsq.sodium: "Na",
        nsq.magnesium: "Mg",
        nsq.aluminum: "Al",
        nsq.silicon: "Si",
        nsq.sulfur: "S",
        nsq.calcium: "Ca",
        nsq.iron: "Fe",
        nsq.nickel: "Ni",
    }
    return element_names.get(code, "Unknown")

def main():
    units = nsq.Const()

    print("=" * 50)
    print("  nuSQuIDS Body Composition Example")
    print("=" * 50)
    print()

    # =========================================================================
    # Part 1: Load Earth model with detailed nuclear composition
    # =========================================================================
    print("--- Part 1: Earth Model with Composition ---")
    print()

    # Load the PREM model with isotopic/nuclear composition
    # The file has 13 columns: r/R, density, ye, H, O, Na, Mg, Al, Si, S, Ca, Fe, Ni
    prem_path = nsq.getResourcePath() + "/astro/EARTH_MODEL_PREM_wIso.dat"
    print(f"Loading Earth model from: {prem_path}")
    print()

    earth = nsq.EarthAtm(prem_path)

    # Create tracks at different zenith angles to probe different depths
    cos_zeniths = [-1.0, -0.8, -0.5, -0.2, 0.0]

    print("Composition at various depths (sampled at track midpoint):")
    print("-" * 80)
    print(f"{'cos(zen)':>10} {'depth(km)':>12} {'density':>10} {'Fe':>8} {'O':>8} {'Si':>8} {'Mg':>8}")
    print("-" * 80)

    for cz in cos_zeniths:
        track = earth.MakeTrackWithCosine(cz)

        # Sample at the midpoint of the track
        x_mid = 0.5 * track.GetFinalX()
        track.SetX(x_mid)

        # Get density and composition at this point
        density = earth.density(track)
        comp = earth.composition(track)

        # Calculate approximate depth
        depth_km = earth.GetRadius() - x_mid / units.km

        # Get element fractions (with default of 0 if not present)
        fe = comp.get(nsq.iron, 0.0)
        o = comp.get(nsq.oxygen, 0.0)
        si = comp.get(nsq.silicon, 0.0)
        mg = comp.get(nsq.magnesium, 0.0)

        print(f"{cz:10.2f} {depth_km:12.2f} {density:10.2f} {fe:8.4f} {o:8.4f} {si:8.4f} {mg:8.4f}")

    print()

    # =========================================================================
    # Part 2: Create ConstantDensity body with custom composition
    # =========================================================================
    print("--- Part 2: ConstantDensity with Custom Composition ---")
    print()

    # Create a water-equivalent material (useful for ice/water detectors)
    water_composition = {
        nsq.hydrogen: 2.0/3.0,  # 2 H atoms per water molecule
        nsq.oxygen: 1.0/3.0,    # 1 O atom per water molecule
    }

    water_density = 1.0  # g/cm^3
    water_ye = 0.555     # electron fraction for water

    water = nsq.ConstantDensity(water_density, water_ye, water_composition)
    water_track = nsq.ConstantDensity.Track(1000.0 * units.km)

    print("Water-equivalent material:")
    print(f"  Density: {water.density(water_track)} g/cm^3")
    print(f"  Electron fraction: {water.ye(water_track)}")
    print("  Composition:")
    water_comp = water.composition(water_track)
    for code, frac in water_comp.items():
        print(f"    {get_element_name(code)}: {frac:.4f}")
    print()

    # =========================================================================
    # Part 3: Run a simple propagation with composition-aware body
    # =========================================================================
    print("--- Part 3: Propagation with Composition ---")
    print()

    # Energy range
    E_min = 1.0 * units.GeV
    E_max = 1.0e6 * units.GeV
    n_energies = 100

    # Create nuSQuIDS instance
    nus = nsq.nuSQUIDS(nsq.logspace(E_min, E_max, n_energies), 3, nsq.NeutrinoType.neutrino, True)

    # Set mixing parameters
    nus.Set_MixingAngle(0, 1, 0.5837)  # theta_12
    nus.Set_MixingAngle(0, 2, 0.1496)  # theta_13
    nus.Set_MixingAngle(1, 2, 0.8552)  # theta_23
    nus.Set_SquareMassDifference(1, 7.42e-5)  # dm21^2
    nus.Set_SquareMassDifference(2, 2.51e-3)  # dm31^2

    # Set initial state (pure muon neutrino)
    inistate = np.zeros((n_energies, 3))
    inistate[:, 1] = 1.0  # Pure muon neutrino
    nus.Set_initial_state(inistate, nsq.Basis.flavor)

    # Use Earth model with composition
    nus.Set_Body(earth)
    nus.Set_Track(earth.MakeTrackWithCosine(-1.0))  # Vertical upgoing

    print("Propagating muon neutrinos through Earth (cos_zen = -1)...")
    nus.EvolveState()
    print("Done!")
    print()

    # Print survival probability at a few energies
    print("Muon neutrino survival probability:")
    print(f"{'Energy (GeV)':>15} {'P(nu_mu->nu_mu)':>15}")

    e_range = nus.GetERange()
    sample_energies = [10.0, 100.0, 1000.0, 10000.0]

    for E in sample_energies:
        # Find closest energy index
        ie = np.argmin(np.abs(e_range / units.GeV - E))
        prob = nus.EvalFlavorAtNode(1, ie, 1)  # nu_mu flavor
        print(f"{E:15.1f} {prob:15.4f}")

    print()
    print("=" * 50)
    print("  Example Complete!")
    print("=" * 50)

if __name__ == "__main__":
    main()
