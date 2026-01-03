# nuSQuIDS Cross Section Data Files

This directory contains neutrino-nucleon and neutrino-nucleus cross section tables
used by nuSQuIDS for calculating neutrino interaction rates during propagation.

## File Format

The HDF5 cross section files follow the format documented in the nuSQuIDS paper
(Appendix D). See the main README for the paper reference.

## Cross Section Families

### CSMS (Cooper-Sarkar, Mertsch, Sarkar)

**Reference:** [arXiv:1106.3723](https://arxiv.org/abs/1106.3723)

These cross sections are based on next-to-leading-order QCD calculations with
HERAPDF1.5 parton distributions.

| File | Description |
|------|-------------|
| `csms.h5` | Combined cross section table |
| `csms_proton.h5` | Proton target cross sections |
| `csms_neutron.h5` | Neutron target cross sections |
| `csms_square.h5` | Isoscalar nucleon cross sections |

**Citation:**
```bibtex
@article{CooperSarkar:2011pa,
    author = "Cooper-Sarkar, Amanda and Mertsch, Philipp and Sarkar, Subir",
    title = "{The high energy neutrino cross-section in the Standard Model and its uncertainty}",
    eprint = "1106.3723",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1007/JHEP08(2011)042",
    journal = "JHEP",
    volume = "08",
    pages = "042",
    year = "2011"
}
```

### NuSigma (from WimpSim)

**Reference:** [WimpSim nucross3 documentation](http://wimpsim.astroparticle.se/code/nucross3.pdf)

These cross sections are computed using the NuSigma library, which is part of the
WimpSim package for WIMP detection simulations.

| File | Description |
|------|-------------|
| `nusigma_sigma_CC.dat` | CC total cross sections (isoscalar) |
| `nusigma_sigma_NC.dat` | NC total cross sections (isoscalar) |
| `nusigma_p_sigma_CC.dat` | CC total cross sections (proton) |
| `nusigma_p_sigma_NC.dat` | NC total cross sections (proton) |
| `nusigma_n_sigma_CC.dat` | CC total cross sections (neutron) |
| `nusigma_n_sigma_NC.dat` | NC total cross sections (neutron) |
| `nusigma_dsde_CC.dat` | CC differential cross sections |
| `nusigma_dsde_NC.dat` | NC differential cross sections |

**Citation:**
```bibtex
@misc{wimpsim,
    author = "Edsjo, Joakim",
    title = "{WimpSim neutrino cross sections}",
    howpublished = "\url{http://wimpsim.astroparticle.se/}",
    note = "See nucross3.pdf documentation"
}
```

### WCG24 (Weigel, Conrad, Garcia 2024)

**Reference:** [arXiv:2408.05866](https://arxiv.org/abs/2408.05866)

These cross sections provide per-nucleus neutrino deep inelastic scattering calculations,
enabling accurate modeling of neutrino propagation through materials with known
nuclear composition. This is particularly important for precision atmospheric
neutrino studies through the Earth.

| File | Target | PDG Code |
|------|--------|----------|
| `wcg24_base_isoscalar.h5` | Isoscalar nucleon | 81 |
| `wcg24_base_proton.h5` | Proton | 2212 |
| `wcg24_base_neutron.h5` | Neutron | 2112 |
| `wcg24_oxygen.h5` | Oxygen-16 | 1000080160 |
| `wcg24_carbon.h5` | Carbon-12 | 1000060120 |
| `wcg24_sodium.h5` | Sodium-23 | 1000110230 |
| `wcg24_magnesium.h5` | Magnesium-24 | 1000120240 |
| `wcg24_aluminum.h5` | Aluminum-27 | 1000130270 |
| `wcg24_silicon.h5` | Silicon-28 | 1000140280 |
| `wcg24_sulfur.h5` | Sulfur-32 | 1000160320 |
| `wcg24_calcium.h5` | Calcium-40 | 1000200400 |
| `wcg24_iron.h5` | Iron-56 | 1000260560 |
| `wcg24_nickel.h5` | Nickel-58 | 1000280580 |
| `wcg24_lead.h5` | Lead-208 | 1000822080 |

**Citation:**
```bibtex
@article{Weigel:2024gzh,
    author = "Weigel, Philip L. R. and Conrad, Janet M. and Garcia-Soto, Alfonso",
    title = "{Cross sections and inelasticity distributions of high-energy neutrino deep inelastic scattering}",
    eprint = "2408.05866",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1103/PhysRevD.111.043044",
    journal = "Phys. Rev. D",
    volume = "111",
    number = "4",
    pages = "043044",
    year = "2025"
}
```

## Usage

### Default Cross Sections

By default, nuSQuIDS loads the CSMS proton/neutron cross sections:

```cpp
// Automatic loading via loadDefaultCrossSections()
nuSQUIDS nus(...);  // Uses csms_proton.h5 and csms_neutron.h5
```

### Loading WCG24 Nuclear Cross Sections

To use per-nucleus cross sections for composition-aware propagation:

```cpp
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/resources.h>

// Create cross section library with nuclear targets
CrossSectionLibrary lib;
std::string xsdir = getResourcePath() + "/xsections/";

// Add per-nucleus cross sections
lib.addTarget(oxygen, NeutrinoDISCrossSectionsFromTables(xsdir + "wcg24_oxygen.h5"));
lib.addTarget(iron, NeutrinoDISCrossSectionsFromTables(xsdir + "wcg24_iron.h5"));
lib.addTarget(silicon, NeutrinoDISCrossSectionsFromTables(xsdir + "wcg24_silicon.h5"));
// ... add other elements as needed

// Add electron target for Glashow resonance
lib.addTarget(electron, GlashowResonanceCrossSection());

// Create nuSQuIDS with custom cross sections
nuSQUIDS nus(...);
nus.Set_CrossSectionLibrary(std::make_shared<CrossSectionLibrary>(lib));
```

### Python Example

```python
import nuSQuIDS as nsq

# Load WCG24 cross sections for specific elements
xsdir = nsq.getResourcePath() + "/xsections/"
lib = nsq.CrossSectionLibrary()
lib.addTarget(nsq.oxygen, nsq.NeutrinoDISCrossSectionsFromTables(xsdir + "wcg24_oxygen.h5"))
lib.addTarget(nsq.iron, nsq.NeutrinoDISCrossSectionsFromTables(xsdir + "wcg24_iron.h5"))
# ... etc

nus = nsq.nuSQUIDS(...)
nus.Set_CrossSectionLibrary(lib)
```

## Acknowledgments

We are grateful to the authors of the cross section calculations for making their
results publicly available. Users of nuSQuIDS are encouraged to cite the appropriate
references when using these cross sections in their work.
