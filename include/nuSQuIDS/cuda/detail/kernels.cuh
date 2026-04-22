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
// Per-thread RK4 correction buffers are sized for MAX_PAIRS (rho, ie)
// pairs; the invariant is nrhos * ceil(ne / EVOLVE_THREADS) <= MAX_PAIRS.
// launchEvolve also checks this at runtime and throws before launch.
// ============================================================
static constexpr int EVOLVE_THREADS = 128;
static constexpr int MAX_PAIRS = 32;

// At nrhos=2 the invariant is 2*ceil(ne/EVOLVE_THREADS) <= MAX_PAIRS,
// i.e. ne <= (MAX_PAIRS/2) * EVOLVE_THREADS = 16 * 128 = 2048.
// Beyond that the runtime check in launchEvolve() rejects the launch.
static_assert(MAX_PAIRS * EVOLVE_THREADS >= 4096,
              "MAX_PAIRS*EVOLVE_THREADS must cover ne*nrhos for the tested grid");

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

/// Main evolution kernel — full propagation per path with adaptive RK4 solver
/// Grid: [1, n_paths] x 128 threads per block
void launchEvolve(const PhysicsParams& params,
                  const PathDeviceData* d_paths,
                  const double* d_H0_array,
                  const double* d_b1_proj,
                  const InteractionDataGPU* d_interaction_data,
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
