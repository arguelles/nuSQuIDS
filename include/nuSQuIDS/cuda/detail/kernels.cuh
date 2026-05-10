/******************************************************************************
*   CUDA kernel declarations for nuSQuIDS CUDA backend.
*   All kernel launch wrappers are declared here.
*   Adapted from CUDAnuSQuIDS (GPL-3.0).
******************************************************************************/

#ifndef NUSQUIDS_CUDA_KERNELS_CUH
#define NUSQUIDS_CUDA_KERNELS_CUH

#include <cuda_runtime.h>
#include "physics.cuh"
#include "body_gpu.cuh"
#include "solver_gpu.cuh"
#include "interactions_gpu.cuh"

namespace nusquids { namespace cuda {

// ============================================================
// Launch-configuration constants shared by the launcher and the kernel.
// RK4 staging buffers live in a per-path slice of a persistent device
// workspace allocated by the Propagator, sized nrhos*ne*SU doubles per path.
// ============================================================
static constexpr int EVOLVE_THREADS = 128;

// ============================================================
// Kernel launch wrappers
// ============================================================

/// Initialize H0 array on device: H0[rho][ie] = dm^2 / (2E)
/// Grid: [ceil(ne/128), 1]
void launchInitH0(const double* d_energies, const double* d_dm2,
                  int ne, int nrhos, int numneu,
                  double* d_H0_array, cudaStream_t stream);

/// Set initial states from flux data
/// Grid: [ceil(ne/128), n_paths]
void launchSetInitialStates(const double* d_flux_data,
                            const double* d_b0_proj,
                            int ne, int nrhos, int numneu,
                            double* d_states,
                            int n_paths, cudaStream_t stream);

/// Main evolution kernel — full propagation per path with adaptive solver.
///
/// Oscillation-only: classic RK4 + step-doubling Richardson extrapolation.
/// Interactions:     Dormand-Prince 5(4) (DOPRI5) embedded pair with FSAL.
///
/// `d_workspace_k` is an array of 7 device pointers, each sized
/// n_paths*nrhos*ne*SU doubles, used to hold k1..k7 stage derivatives for
/// the interaction-path DOPRI5 integrator. May be nullptr when interactions
/// are disabled (the kernel will not access them in that case).
///
/// Grid: [1, n_paths] x 128 threads per block
void launchEvolve(const PhysicsParams& params,
                  const PathDeviceData* d_paths,
                  const double* d_H0_array,
                  const double* d_b1_proj,
                  const InteractionDataGPU* d_interaction_data,
                  double* d_workspace_corrected,
                  double* d_workspace_sf,
                  double* const* d_workspace_k,   ///< [7] k-stage slabs (interaction path)
                  const SolverConfig& solver_config,
                  double* d_states,
                  int n_paths, int numneu,
                  cudaStream_t stream);

/// Extract flavor expectation values from final states
/// Grid: [ceil(ne/128), n_paths]
void launchEvalFlavors(const double* d_states,
                       const double* d_evol_proj,
                       int ne, int nrhos, int numneu,
                       double* d_flavors,
                       int n_paths, cudaStream_t stream);

}} // namespace nusquids::cuda

#endif // NUSQUIDS_CUDA_KERNELS_CUH
