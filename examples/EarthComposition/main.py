#!/usr/bin/env python3
"""
EarthComposition Example - Comparing WCG24 Proton/Neutron vs Nuclear Cross Sections

This example demonstrates the difference between:
1. WCG24 proton/neutron mode: Per-nucleon cross sections (wcg24_base_proton/neutron.h5)
2. WCG24 nuclear mode: Per-element cross sections with PREM composition

Both use the same WCG24 cross section family (arXiv:2408.05866) for a fair comparison.
The nuclear mode uses detailed Earth composition from EARTH_MODEL_PREM_wIso.dat.

We propagate atmospheric neutrinos through Earth using nuSQUIDSAtm and compare
the final flux ratios for both modes.
"""

import numpy as np
import matplotlib.pyplot as plt
import nuSQuIDS as nsq

# Units
GeV = nsq.Const().GeV
km = nsq.Const().km


def run_propagation(use_nuclear_xs=False):
    """
    Run atmospheric neutrino propagation.

    Parameters
    ----------
    use_nuclear_xs : bool
        If True, use WCG24 per-nucleus cross sections with PREM composition.
        If False, use WCG24 proton/neutron cross sections with isoscalar PREM.

    Returns
    -------
    tuple
        (energies, cos_zeniths, flux_ratios) where flux_ratios[icz, ie, flv]
        is the ratio of final to initial flux for each zenith, energy, and flavor.
    """
    # Energy range: 10 GeV to 1 PeV
    Emin = 10.0 * GeV
    Emax = 1.0e6 * GeV
    n_energies = 100

    # Zenith angles: upgoing through horizontal
    cos_zenith_min = -1.0
    cos_zenith_max = 0.0
    n_cos_zenith = 20

    cos_zeniths = nsq.linspace(cos_zenith_min, cos_zenith_max, n_cos_zenith)
    energies = nsq.logspace(Emin, Emax, n_energies)

    # Create nuSQUIDSAtm with interactions
    nus_atm = nsq.nuSQUIDSAtm(cos_zeniths, energies, 3, nsq.NeutrinoType.both, True)

    # Set cross sections and Earth model
    if use_nuclear_xs:
        # Load WCG24 per-nucleus cross sections
        xslib = nsq.loadWCG24NuclearCrossSections()
        nus_atm.SetNeutrinoCrossSections(xslib)

        # Load PREM with nuclear composition (13-column file)
        prem_path = nsq.getResourcePath() + "/astro/EARTH_MODEL_PREM_wIso.dat"
        earth = nsq.EarthAtm(prem_path)
        mode_name = "nuclear (WCG24 per-element)"
    else:
        # Load WCG24 proton/neutron cross sections
        xslib = nsq.loadWCG24CrossSections()
        nus_atm.SetNeutrinoCrossSections(xslib)

        # Use default isoscalar PREM (3-column file)
        earth = nsq.EarthAtm()
        mode_name = "nucleon (WCG24 p/n)"

    nus_atm.Set_EarthModel(earth)

    # Set mixing parameters (NuFit 5.0 values)
    nus_atm.Set_MixingAngle(0, 1, 0.5837)   # theta_12
    nus_atm.Set_MixingAngle(0, 2, 0.1496)   # theta_13
    nus_atm.Set_MixingAngle(1, 2, 0.8552)   # theta_23
    nus_atm.Set_SquareMassDifference(1, 7.42e-5)   # dm21^2
    nus_atm.Set_SquareMassDifference(2, 2.51e-3)   # dm31^2
    nus_atm.Set_CPPhase(0, 2, 0.0)

    # Set numerical parameters
    nus_atm.Set_rel_error(1.0e-6)
    nus_atm.Set_abs_error(1.0e-6)
    nus_atm.Set_GSL_step(nsq.GSL_STEP_FUNCTIONS.GSL_STEP_RK4)

    # Set initial state: E^-2 spectrum for all flavors (atmospheric-like)
    # Initial state shape: [n_cos_zenith, n_energies, 2 (nu/nubar), 3 (flavors)]
    initial_state = np.zeros((n_cos_zenith, n_energies, 2, 3))

    E_array = np.array([energies[i] / GeV for i in range(n_energies)])

    for icz in range(n_cos_zenith):
        for ie in range(n_energies):
            E = E_array[ie]
            flux = E**(-2)  # E^-2 spectrum
            for rho in range(2):  # neutrino and antineutrino
                # Atmospheric ratio approximately: nu_e : nu_mu : nu_tau = 1 : 2 : 0
                initial_state[icz, ie, rho, 0] = flux * 1.0  # nu_e
                initial_state[icz, ie, rho, 1] = flux * 2.0  # nu_mu
                initial_state[icz, ie, rho, 2] = flux * 0.0  # nu_tau

    nus_atm.Set_initial_state(initial_state, nsq.Basis.flavor)
    nus_atm.Set_IncludeOscillations(True)

    # Evolve
    print(f"Propagating {mode_name}...")
    nus_atm.EvolveState()
    print("Done!")

    # Extract final fluxes and compute ratios
    flux_ratios = np.zeros((n_cos_zenith, n_energies, 3))  # [cz, E, flavor]

    for icz in range(n_cos_zenith):
        for ie in range(n_energies):
            for flv in range(3):
                # Sum over neutrino and antineutrino
                final_flux = (nus_atm.EvalFlavor(flv, cos_zeniths[icz], energies[ie], 0) +
                              nus_atm.EvalFlavor(flv, cos_zeniths[icz], energies[ie], 1))
                initial_flux = initial_state[icz, ie, 0, flv] + initial_state[icz, ie, 1, flv]

                if initial_flux > 0:
                    flux_ratios[icz, ie, flv] = final_flux / initial_flux
                else:
                    # For tau (initial flux = 0), store relative appearance
                    # normalized to initial mu flux for comparison
                    mu_initial = initial_state[icz, ie, 0, 1] + initial_state[icz, ie, 1, 1]
                    if mu_initial > 0:
                        flux_ratios[icz, ie, flv] = final_flux / mu_initial
                    else:
                        flux_ratios[icz, ie, flv] = 0.0

    return E_array, np.array([cos_zeniths[i] for i in range(n_cos_zenith)]), flux_ratios


def main():
    print("=" * 70)
    print("EarthComposition Example")
    print("Comparing WCG24 proton/neutron vs per-element cross sections")
    print("Reference: arXiv:2408.05866")
    print("=" * 70)
    print()

    # Run both propagations
    E_pn, cz_pn, ratios_pn = run_propagation(use_nuclear_xs=False)
    E_nuc, cz_nuc, ratios_nuc = run_propagation(use_nuclear_xs=True)

    # Create figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    flavor_names = [r'$\nu_e$', r'$\nu_\mu$', r'$\nu_\tau$']

    # Select a few zenith angles to plot
    cz_indices = [0, 5, 10, 15, 19]  # -1.0, -0.75, -0.5, -0.25, 0.0 approximately

    # Top row: Flux ratios for each flavor (p/n vs nuclear)
    # Note: For tau, initial flux is 0, so we show absolute appearance flux instead
    for flv in range(3):
        ax = axes[0, flv]

        for idx in cz_indices:
            cz = cz_pn[idx]
            label = f'cos(z)={cz:.2f}'

            # Proton/neutron (solid lines)
            ax.semilogx(E_pn, ratios_pn[idx, :, flv], '-',
                       label=f'{label} (p/n)', alpha=0.7)
            # Nuclear (dashed lines)
            ax.semilogx(E_nuc, ratios_nuc[idx, :, flv], '--',
                       label=f'{label} (nuc)', alpha=0.7)

        ax.set_xlabel('Energy [GeV]')
        if flv == 2:  # tau
            ax.set_ylabel(r'$\nu_\tau^{final} / \nu_\mu^{initial}$')
            ax.set_title(r'$\nu_\tau$ Appearance')
        else:
            ax.set_ylabel(f'{flavor_names[flv]} Flux Ratio (Final/Initial)')
            ax.set_title(f'{flavor_names[flv]} Survival/Appearance')
        ax.set_xlim(10, 1e6)
        ax.set_ylim(0, 1.5)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=6, ncol=2)

    # Bottom row: Ratio of nuclear to p/n for each flavor
    for flv in range(3):
        ax = axes[1, flv]

        for idx in cz_indices:
            cz = cz_pn[idx]

            # Avoid division by zero
            ratio = np.divide(ratios_nuc[idx, :, flv], ratios_pn[idx, :, flv],
                            out=np.ones_like(ratios_pn[idx, :, flv]),
                            where=ratios_pn[idx, :, flv] > 1e-10)

            ax.semilogx(E_pn, ratio, '-', label=f'cos(z)={cz:.2f}')

        ax.set_xlabel('Energy [GeV]')
        ax.set_ylabel('Nuclear / (p/n)')
        ax.set_title(f'{flavor_names[flv]} Nuclear Composition Effect')
        ax.set_xlim(10, 1e6)
        ax.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
        ax.set_ylim(0.8, 1.2)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    plt.suptitle('Effect of Nuclear Composition on Atmospheric Neutrino Propagation\n'
                 'WCG24 cross sections (arXiv:2408.05866), Initial spectrum: $E^{-2}$', fontsize=14)
    plt.tight_layout()
    plt.savefig('earth_composition_comparison.png', dpi=150, bbox_inches='tight')
    plt.savefig('earth_composition_comparison.eps', bbox_inches='tight')
    print()
    print("Saved: earth_composition_comparison.png and earth_composition_comparison.eps")
    plt.show()


if __name__ == "__main__":
    main()
