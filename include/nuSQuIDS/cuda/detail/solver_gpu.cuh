/******************************************************************************
*   GPU-resident ODE solver for nuSQuIDS CUDA backend.
*   RK4 with step-doubling and adaptive step size control.
*   Adapted from CUDAnuSQuIDS Version2 solver (GPL-3.0).
*
*   One thread block per path, threads parallelize over energies.
*   Uses warp-level reductions for error estimation (sm_80+).
******************************************************************************/

#ifndef NUSQUIDS_CUDA_SOLVER_GPU_CUH
#define NUSQUIDS_CUDA_SOLVER_GPU_CUH

#include <cuda_runtime.h>
#include "physics.cuh"
#include "body_gpu.cuh"

namespace nusquids { namespace cuda {

// ============================================================
// Solver configuration
// ============================================================

struct SolverConfig {
  double h_min;          ///< Minimum step size
  double h_max;          ///< Maximum step size
  double h_initial;      ///< Initial step size
  double rel_error;      ///< Relative error tolerance
  double abs_error;       ///< Absolute error tolerance
  int max_steps;         ///< Maximum number of steps per path
};

// ============================================================
// RK4 step functions
// ============================================================

/// Perform one classical RK4 step for all energy nodes handled by this thread
/// \param state Current state array [ne * nrhos * su_size]
/// \param h Step size
/// \param x Current position
/// \param result Output state after step
/// All parameters reference device memory
template<int NFLV>
__device__
void rk4Step(const double* __restrict__ state,
             double h, double x,
             const PhysicsParams& params,
             const PathDeviceData& path,
             const double* __restrict__ H0_array,
             const double* __restrict__ b1_proj,
             const double* __restrict__ interaction_data,
             double* __restrict__ k_buffer,
             double* __restrict__ result);

/// Compute the RHS of the evolution equation for one energy node
/// dρ/dt = -i[H, ρ] + D[ρ]
/// where D[ρ] includes absorption and cascading interactions
template<int NFLV>
__device__
void computeRHS(int ie, int rho,
                double x, double density, double ye,
                const double* __restrict__ state,
                const PhysicsParams& params,
                const double* __restrict__ H0_array,
                const double* __restrict__ evol_proj,
                const double* __restrict__ interaction_data,
                double* __restrict__ result);

// ============================================================
// Adaptive step size control
// ============================================================

/// Compute error estimate from step-doubling:
/// Do one full step of size h and two half steps of size h/2,
/// compare results to estimate error.
/// Returns the maximum relative error across all energy nodes in the block.
template<int NFLV>
__device__
double estimateError(const double* __restrict__ full_step,
                     const double* __restrict__ half_step,
                     int ne, int nrhos, int su_size);

/// Block-wide maximum reduction using warp shuffles (sm_80+)
__device__ __forceinline__
double blockReduceMax(double val) {
  // Warp-level reduction
  for (int offset = warpSize / 2; offset > 0; offset /= 2)
    val = fmax(val, __shfl_down_sync(0xffffffff, val, offset));

  // Inter-warp reduction via shared memory
  __shared__ double warp_max[32]; // max 32 warps per block
  int lane = threadIdx.x % warpSize;
  int warp_id = threadIdx.x / warpSize;

  if (lane == 0)
    warp_max[warp_id] = val;
  __syncthreads();

  // First warp reduces across warps
  int num_warps = (blockDim.x + warpSize - 1) / warpSize;
  val = (threadIdx.x < num_warps) ? warp_max[threadIdx.x] : 0.0;
  if (warp_id == 0) {
    for (int offset = warpSize / 2; offset > 0; offset /= 2)
      val = fmax(val, __shfl_down_sync(0xffffffff, val, offset));
  }

  __shared__ double block_max;
  if (threadIdx.x == 0)
    block_max = val;
  __syncthreads();
  return block_max;
}

// ============================================================
// Main evolution kernel
// ============================================================

/// Main monolithic evolution kernel.
/// Each block handles one path (zenith angle).
/// Threads within a block parallelize over energy nodes.
///
/// The kernel:
/// 1. Loads path's density profile spline data
/// 2. Steps through time with adaptive RK4 (step-doubling)
/// 3. At each step: updates track position → evaluates density spline →
///    evolves projectors → computes interactions → evaluates RHS → advances state
/// 4. Writes final state to global memory
///
/// Launch: <<<n_paths, 128>>>
template<int NFLV>
__global__
void __launch_bounds__(128, 4)
evolveKernel(const PhysicsParams params,
             const PathDeviceData* __restrict__ paths,
             const double* __restrict__ H0_array,
             const double* __restrict__ b1_proj,
             const double* __restrict__ interaction_data,
             const SolverConfig solver_config,
             double* __restrict__ states,
             int n_paths);

/// Kernel to initialize H0 array: H0[rho][ie] = dm^2 / (2E)
template<int NFLV>
__global__
void initH0Kernel(const double* __restrict__ energies,
                  const double* __restrict__ dm2,
                  int ne, int nrhos,
                  double* __restrict__ H0_array);

/// Kernel to set initial states from uploaded flux data
template<int NFLV>
__global__
void setInitialStatesKernel(const double* __restrict__ flux_data,
                            const double* __restrict__ b0_proj,
                            int ne, int nrhos, int numneu,
                            double* __restrict__ states,
                            int n_paths);

}} // namespace nusquids::cuda

#endif // NUSQUIDS_CUDA_SOLVER_GPU_CUH
