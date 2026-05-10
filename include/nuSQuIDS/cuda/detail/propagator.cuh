/******************************************************************************
*   Propagator class for nuSQuIDS CUDA backend.
*   Manages one GPU device: memory allocation, data transfer, batch evolution.
*   Adapted from CUDAnuSQuIDS (GPL-3.0).
******************************************************************************/

#ifndef NUSQUIDS_CUDA_PROPAGATOR_CUH
#define NUSQUIDS_CUDA_PROPAGATOR_CUH

#include <cuda_runtime.h>
#include <vector>
#include <memory>

#include "memory.cuh"
#include "physics.cuh"
#include "body_gpu.cuh"
#include "interactions_gpu.cuh"
#include "solver_gpu.cuh"

namespace nusquids { namespace cuda {

/// Manages one GPU device for neutrino propagation
class Propagator {
public:
  /// \param device_id CUDA device to use
  /// \param batch_size_limit Maximum paths per batch
  Propagator(int device_id, int batch_size_limit);
  ~Propagator();

  // Non-copyable
  Propagator(const Propagator&) = delete;
  Propagator& operator=(const Propagator&) = delete;

  /// Upload shared physics data (done once per evolution)
  /// This includes H0 array, projectors, interaction data
  /// \param interaction_host_data  Pointer to InteractionDataHost (from cuda_backend.h),
  ///                               or nullptr for oscillation-only. Passed as void* to
  ///                               avoid coupling this header to cuda_backend.h.
  void uploadSharedData(const PhysicsParams& params,
                        const void* interaction_host_data,
                        const double* H0_array_host,
                        const double* b1_proj_host,
                        const SolverConfig& solver_config);

  /// Evolve a batch of paths
  /// \param profiles CPU-side density profiles for each path
  /// \param initial_states Initial state data [n_paths * ne * nrhos * su_size]
  /// \param final_states Output: evolved states [same layout]
  /// \param n_paths Number of paths in this batch
  void evolveBatch(const std::vector<GPUDensityProfile>& profiles,
                   const double* initial_states,
                   double* final_states,
                   int n_paths);

  int getDeviceId() const { return device_id_; }

private:
  int device_id_;
  int batch_size_limit_;

  // CUDA resources
  cudaStream_t stream_;
  cudaEvent_t event_;

  // Device memory for current batch
  double* d_states_;           ///< State data on device
  PathDeviceData* d_paths_;    ///< Per-path data on device

  // Persistent integrator staging workspaces (grow-on-demand, freed in destructor).
  // Sized n_paths * nrhos * ne * SU doubles each; replace the per-thread
  // corrected_buf/sf_buf locals that dominated stack spill.
  //
  // For oscillation-only (RK4 + Richardson) only `corrected` and `sf` are used.
  // For the interaction path we use Dormand-Prince 5(4) and need 7 additional
  // slabs to hold k1..k7. They are allocated lazily the first time an
  // interaction batch is evolved.
  double* d_workspace_corrected_;
  double* d_workspace_sf_;
  double* d_workspace_k_[7];     ///< k1..k7 stage derivatives (DOPRI5)
  size_t workspace_size_bytes_;
  size_t workspace_k_size_bytes_; ///< current size of each k slab

  // Shared data (persists across batches)
  double* d_H0_array_;
  double* d_b1_proj_;
  PhysicsParams params_;
  SolverConfig solver_config_;
  InteractionDataGPU interaction_data_;

  // Spline data buffer for profiles
  double* d_spline_buffer_;
  size_t spline_buffer_size_;

  bool shared_data_uploaded_;

  /// Upload profile spline data for a batch of paths
  void uploadProfiles(const std::vector<GPUDensityProfile>& profiles,
                      std::vector<PathDeviceData>& host_paths);

  /// Select the appropriate templated kernel based on numneu
  void launchEvolveKernel(int n_paths);

  /// Allocate device memory for batch
  void allocateBatch(int n_paths);
  void freeBatch();
};

}} // namespace nusquids::cuda

#endif // NUSQUIDS_CUDA_PROPAGATOR_CUH
