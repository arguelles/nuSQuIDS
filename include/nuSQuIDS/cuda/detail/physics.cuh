/******************************************************************************
*   Device-side physics operators for nuSQuIDS CUDA backend.
*   H0, HI, GammaRho, InteractionsRho on GPU.
*   Adapted from CUDAnuSQuIDS (GPL-3.0).
******************************************************************************/

#ifndef NUSQUIDS_CUDA_PHYSICS_CUH
#define NUSQUIDS_CUDA_PHYSICS_CUH

#include <cuda_runtime.h>
#include "body_gpu.cuh"

namespace nusquids { namespace cuda {

// ============================================================
// Physical constants (natural units, matching SQuIDS)
// ============================================================

// Fermi constant times Avogadro's number [eV^-2 * cm^3 * g^-1]
// GF * Na / sqrt(2), from SQuIDS Const
__device__ constexpr double SQRT2_GF_NA = 7.63247e-14; // sqrt(2) * GF * Na

// ============================================================
// Device-side physics parameters (uploaded once)
// ============================================================

/// Physics parameters shared across all paths
struct PhysicsParams {
  int numneu;           ///< Number of neutrino flavors
  int nrhos;            ///< Number of density matrix equations (1 or 2)
  int ne;               ///< Number of energy nodes
  int su_size;          ///< Size of SU(N) vector: numneu^2

  // Mixing parameters (used to build H0)
  // dm2[i] = mass-squared difference dm^2_{i1} for i=1..numneu-1
  double dm2[6];        ///< Mass squared differences (max 6 for numneu=4)
  int n_dm2;

  // H0 is pre-computed on CPU and uploaded per energy node
  // H0[rho][ie] = numneu^2 components

  // Matter potential constant: sqrt(2) * GF * Na * cm^{-3} [eV]
  double HI_constants;

  // Neutrino type: 1=neutrino, 2=antineutrino, 3=both
  int NT_type;

  // Interaction flags
  bool iinteraction;
  bool ioscillations;
  bool tauregeneration;
  bool iglashow;

  // Basis
  int basis;            ///< 0=mass, 1=flavor, 2=interaction

  // Unit conversion for interaction densities
  // density_natural = density_g_cm3 * gr_to_eV_cm3
  // number_density = density_natural / nucleon_mass
  double gr_to_eV_cm3;  ///< params.gr * pow(params.cm, -3) [eV^4 per g/cm³]
  double proton_mass;    ///< [eV]
  double neutron_mass;   ///< [eV]
  double electron_mass;  ///< [eV]
};

/// Per-path device data: state + evolved projectors + body profile
struct PathDeviceData {
  GPUDensityProfileDevice profile;
  double xini;          ///< Track initial position
  double xend;          ///< Track final position
  double time_offset;   ///< Time offset between SQuIDS time and Track(x)
};

// ============================================================
// Device-side H0 computation
// ============================================================

/// Compute H0 for a single energy node
/// H0(E, rho) = dm2 / (2E) in the appropriate basis
/// \param dm2 Mass squared differences
/// \param E Energy in eV
/// \param rho 0=neutrino, 1=antineutrino
/// \param numneu Number of flavors
/// \param result Output SU(N) vector
template<int NFLV>
__device__
void computeH0(const double* __restrict__ dm2, double E, int rho,
               double* __restrict__ result);

// ============================================================
// Device-side matter potential (HI)
// ============================================================

/// Compute the matter potential at current position
/// V_CC = sqrt(2) * GF * Ne for electron neutrinos
/// V_NC = -GF * Nn / sqrt(2) for all flavors
/// \param density Matter density [g/cm^3]
/// \param ye Electron fraction
/// \param numneu Number of flavors
/// \param rho 0=neutrino, 1=antineutrino
/// \param evol_proj Evolved flavor projectors at current time [nflv][su_size]
/// \param result Output SU(N) vector for matter potential
template<int NFLV>
__device__
void computeHI(double density, double ye, int rho,
               const double* __restrict__ evol_proj,
               double* __restrict__ result);

// ============================================================
// Device-side absorption (GammaRho)
// ============================================================

/// Compute absorption term from inverse interaction lengths
/// GammaRho = -0.5 * sum_flv invlen[flv] * evol_proj[flv]
template<int NFLV>
__device__
void computeGammaRho(int ie, int rho, int ne,
                     const double* __restrict__ invlen_CC,
                     const double* __restrict__ invlen_NC,
                     const double* __restrict__ invlen_GR,
                     const double* __restrict__ evol_proj,
                     double* __restrict__ result);

// ============================================================
// Device-side cascading interactions (InteractionsRho)
// ============================================================

/// Compute NC cascading + tau regeneration + Glashow resonance
/// This is the O(ne^2) kernel that computes energy redistribution
template<int NFLV>
__device__
void computeInteractionsRho(int ie, int rho, int ne,
                            const double* __restrict__ state,
                            const double* __restrict__ dNdE_NC,
                            const double* __restrict__ dNdE_CC,
                            const double* __restrict__ dNdE_tau,
                            const double* __restrict__ dNdE_GR,
                            const double* __restrict__ delE,
                            const double* __restrict__ evol_proj,
                            const double* __restrict__ target_fractions,
                            int n_targets,
                            bool tauregeneration,
                            bool iglashow,
                            double* __restrict__ result);

// ============================================================
// Projector evolution
// ============================================================

/// Evolve flavor projectors under H0 to time t
/// evol_proj[flv][ie] = exp(-iH0*t) * proj[flv] * exp(iH0*t)
template<int NFLV>
__device__
void evolveProjectors(double time, int ie, int rho,
                      const double* __restrict__ H0,
                      const double* __restrict__ b1_proj,
                      double* __restrict__ evol_proj);

// ============================================================
// Interaction state update
// ============================================================

/// Update inverse interaction lengths from cross sections × density
/// invlen = density * Na * sigma * target_fraction
template<int NFLV>
__device__
void updateInteractionLengths(double density, double ye,
                              int ie, int rho, int ne,
                              const double* __restrict__ sigma_CC,
                              const double* __restrict__ sigma_NC,
                              const double* __restrict__ sigma_GR,
                              const double* __restrict__ target_fractions,
                              int n_targets,
                              double* __restrict__ invlen_CC,
                              double* __restrict__ invlen_NC,
                              double* __restrict__ invlen_GR);

}} // namespace nusquids::cuda

#endif // NUSQUIDS_CUDA_PHYSICS_CUH
