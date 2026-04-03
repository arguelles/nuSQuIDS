/******************************************************************************
*   Propagator implementation for nuSQuIDS CUDA backend.
*   Manages one GPU device: memory allocation, data transfer, batch evolution.
*   Adapted from CUDAnuSQuIDS (GPL-3.0).
******************************************************************************/

#include "nuSQuIDS/cuda/detail/propagator.cuh"
#include "nuSQuIDS/cuda/detail/kernels.cuh"
#include "nuSQuIDS/cuda/detail/memory.cuh"
#include "nuSQuIDS/cuda/detail/interactions_gpu.cuh"
#include "nuSQuIDS/cuda/cuda_backend.h"  // InteractionDataHost

#include <algorithm>
#include <cstring>

namespace nusquids { namespace cuda {

// ---------------------------------------------------------------------------
// InteractionDataGPU management (declared in interactions_gpu.cuh)
// ---------------------------------------------------------------------------

void freeInteractionData(InteractionDataGPU& data) {
  if (data.d_dNdE_CC)      { cudaFree(data.d_dNdE_CC);      data.d_dNdE_CC = nullptr; }
  if (data.d_dNdE_NC)      { cudaFree(data.d_dNdE_NC);      data.d_dNdE_NC = nullptr; }
  if (data.d_sigma_CC)     { cudaFree(data.d_sigma_CC);     data.d_sigma_CC = nullptr; }
  if (data.d_sigma_NC)     { cudaFree(data.d_sigma_NC);     data.d_sigma_NC = nullptr; }
  if (data.d_dNdE_tau_all) { cudaFree(data.d_dNdE_tau_all); data.d_dNdE_tau_all = nullptr; }
  if (data.d_dNdE_tau_lep) { cudaFree(data.d_dNdE_tau_lep); data.d_dNdE_tau_lep = nullptr; }
  if (data.d_sigma_GR)     { cudaFree(data.d_sigma_GR);     data.d_sigma_GR = nullptr; }
  if (data.d_dNdE_GR)      { cudaFree(data.d_dNdE_GR);      data.d_dNdE_GR = nullptr; }
  if (data.d_energies)     { cudaFree(data.d_energies);     data.d_energies = nullptr; }
  if (data.d_delE)         { cudaFree(data.d_delE);         data.d_delE = nullptr; }
  if (data.d_b0_proj)      { cudaFree(data.d_b0_proj);      data.d_b0_proj = nullptr; }
  if (data.d_b1_proj)      { cudaFree(data.d_b1_proj);      data.d_b1_proj = nullptr; }
  data.total_bytes = 0;
}

InteractionDataGPU uploadInteractionData(const void* int_data_host_ptr,
                                         int ne, int nrhos, int numneu,
                                         const double* energies,
                                         const double* delE,
                                         const double* b0_proj,
                                         const double* b1_proj,
                                         cudaStream_t stream) {
  InteractionDataGPU data;
  memset(&data, 0, sizeof(data));
  data.ne = ne;
  data.nrhos = nrhos;
  data.numneu = numneu;
  data.n_targets = 0;
  data.total_bytes = 0;

  // If interaction data provides energies/delE, use those as fallbacks
  const nusquids::InteractionDataHost* idata_ptr = nullptr;
  if (int_data_host_ptr)
    idata_ptr = reinterpret_cast<const nusquids::InteractionDataHost*>(int_data_host_ptr);
  const double* ene_src = energies ? energies : (idata_ptr ? idata_ptr->energies : nullptr);
  const double* delE_src = delE ? delE : (idata_ptr ? idata_ptr->delE : nullptr);

  // Upload energy grid
  if (ene_src) {
    size_t ene_bytes = ne * sizeof(double);
    NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_energies, ene_bytes));
    NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_energies, ene_src, ene_bytes,
                                         cudaMemcpyHostToDevice, stream));
    data.total_bytes += ene_bytes;
  }

  if (delE_src && ne > 1) {
    size_t del_bytes = (ne - 1) * sizeof(double);
    NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_delE, del_bytes));
    NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_delE, delE_src, del_bytes,
                                         cudaMemcpyHostToDevice, stream));
    data.total_bytes += del_bytes;
  }

  // Upload projectors
  int su_size = numneu * numneu;
  if (b0_proj) {
    size_t b0_bytes = numneu * su_size * sizeof(double);
    NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_b0_proj, b0_bytes));
    NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_b0_proj, b0_proj, b0_bytes,
                                         cudaMemcpyHostToDevice, stream));
    data.total_bytes += b0_bytes;
  }
  if (b1_proj) {
    size_t b1_bytes = nrhos * numneu * su_size * sizeof(double);
    NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_b1_proj, b1_bytes));
    NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_b1_proj, b1_proj, b1_bytes,
                                         cudaMemcpyHostToDevice, stream));
    data.total_bytes += b1_bytes;
  }

  // Upload cross-section tables if interaction data is provided
  if (idata_ptr) {
    const auto& idata = *idata_ptr;

    data.n_targets = idata.n_targets;

    if (idata.n_targets > 0 && idata.sigma_CC && idata.sigma_NC) {
      // Total cross sections: [n_targets * nrhos * numneu * ne]
      size_t sigma_count = (size_t)idata.n_targets * nrhos * numneu * ne;
      size_t sigma_bytes = sigma_count * sizeof(double);
      NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_sigma_CC, sigma_bytes));
      NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_sigma_CC, idata.sigma_CC, sigma_bytes,
                                           cudaMemcpyHostToDevice, stream));
      data.total_bytes += sigma_bytes;

      NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_sigma_NC, sigma_bytes));
      NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_sigma_NC, idata.sigma_NC, sigma_bytes,
                                           cudaMemcpyHostToDevice, stream));
      data.total_bytes += sigma_bytes;

      // Differential cross sections: [n_targets * nrhos * numneu * ne * ne]
      size_t dNdE_count = sigma_count * ne;
      size_t dNdE_bytes = dNdE_count * sizeof(double);
      if (idata.dNdE_CC) {
        NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_dNdE_CC, dNdE_bytes));
        NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_dNdE_CC, idata.dNdE_CC, dNdE_bytes,
                                             cudaMemcpyHostToDevice, stream));
        data.total_bytes += dNdE_bytes;
      }
      if (idata.dNdE_NC) {
        NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_dNdE_NC, dNdE_bytes));
        NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_dNdE_NC, idata.dNdE_NC, dNdE_bytes,
                                             cudaMemcpyHostToDevice, stream));
        data.total_bytes += dNdE_bytes;
      }
    }

    // Glashow resonance
    if (idata.has_glashow && idata.sigma_GR && idata.dNdE_GR) {
      size_t gr_sigma_bytes = ne * sizeof(double);
      NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_sigma_GR, gr_sigma_bytes));
      NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_sigma_GR, idata.sigma_GR, gr_sigma_bytes,
                                           cudaMemcpyHostToDevice, stream));
      data.total_bytes += gr_sigma_bytes;

      size_t gr_dNdE_bytes = (size_t)ne * ne * sizeof(double);
      NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_dNdE_GR, gr_dNdE_bytes));
      NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_dNdE_GR, idata.dNdE_GR, gr_dNdE_bytes,
                                           cudaMemcpyHostToDevice, stream));
      data.total_bytes += gr_dNdE_bytes;
    }

    // Tau decay spectra
    if (idata.has_tau_regen && idata.dNdE_tau_all && idata.dNdE_tau_lep) {
      size_t tau_bytes = (size_t)nrhos * ne * ne * sizeof(double);
      NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_dNdE_tau_all, tau_bytes));
      NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_dNdE_tau_all, idata.dNdE_tau_all, tau_bytes,
                                           cudaMemcpyHostToDevice, stream));
      data.total_bytes += tau_bytes;

      NUSQUIDS_CUDA_CHECK(cudaMalloc(&data.d_dNdE_tau_lep, tau_bytes));
      NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(data.d_dNdE_tau_lep, idata.dNdE_tau_lep, tau_bytes,
                                           cudaMemcpyHostToDevice, stream));
      data.total_bytes += tau_bytes;
    }
  }

  return data;
}

Propagator::Propagator(int device_id, int batch_size_limit)
  : device_id_(device_id),
    batch_size_limit_(batch_size_limit),
    d_states_(nullptr),
    d_paths_(nullptr),
    d_H0_array_(nullptr),
    d_b1_proj_(nullptr),
    d_spline_buffer_(nullptr),
    spline_buffer_size_(0),
    shared_data_uploaded_(false) {

  NUSQUIDS_CUDA_CHECK(cudaSetDevice(device_id_));

  // The interacting kernel uses deep call stacks (anticommutator, cascade, etc.)
  // Default thread stack is 1024 bytes; we need ~4KB for the full interaction chain
  // Increase stack size for the interaction kernel's deep call chain.
  // This may fail on MIG partitions — that's OK if the kernel doesn't
  // actually need the extra stack (oscillation-only uses less).
  cudaDeviceSetLimit(cudaLimitStackSize, 8192);
  // Always clear any error from the above so it doesn't poison kernel launches
  cudaGetLastError();

  NUSQUIDS_CUDA_CHECK(cudaStreamCreate(&stream_));
  NUSQUIDS_CUDA_CHECK(cudaEventCreate(&event_));
}

Propagator::~Propagator() {
  cudaSetDevice(device_id_);
  freeBatch();
  if (d_H0_array_) cudaFree(d_H0_array_);
  if (d_b1_proj_) cudaFree(d_b1_proj_);
  if (d_spline_buffer_) cudaFree(d_spline_buffer_);
  freeInteractionData(interaction_data_);
  cudaStreamDestroy(stream_);
  cudaEventDestroy(event_);
}

void Propagator::uploadSharedData(const PhysicsParams& params,
                                  const void* interaction_host_data,
                                  const double* H0_array_host,
                                  const double* b1_proj_host,
                                  const SolverConfig& solver_config) {
  NUSQUIDS_CUDA_CHECK(cudaSetDevice(device_id_));

  params_ = params;
  solver_config_ = solver_config;

  // Upload interaction data to this GPU if provided
  freeInteractionData(interaction_data_);
  if (interaction_host_data) {
    interaction_data_ = uploadInteractionData(
      interaction_host_data, params.ne, params.nrhos, params.numneu,
      nullptr, nullptr, nullptr, nullptr, stream_);
  } else {
    memset(&interaction_data_, 0, sizeof(interaction_data_));
  }

  int su_size = params.numneu * params.numneu;

  // Upload H0 array
  size_t H0_bytes = params.nrhos * params.ne * su_size * sizeof(double);
  if (d_H0_array_) cudaFree(d_H0_array_);
  NUSQUIDS_CUDA_CHECK(cudaMalloc(&d_H0_array_, H0_bytes));
  NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(d_H0_array_, H0_array_host, H0_bytes,
                                       cudaMemcpyHostToDevice, stream_));

  // Upload b1 projectors
  size_t proj_bytes = params.nrhos * params.numneu * su_size * sizeof(double);
  if (d_b1_proj_) cudaFree(d_b1_proj_);
  NUSQUIDS_CUDA_CHECK(cudaMalloc(&d_b1_proj_, proj_bytes));
  NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(d_b1_proj_, b1_proj_host, proj_bytes,
                                       cudaMemcpyHostToDevice, stream_));

  shared_data_uploaded_ = true;
}

void Propagator::allocateBatch(int n_paths) {
  NUSQUIDS_CUDA_CHECK(cudaSetDevice(device_id_));

  int su_size = params_.numneu * params_.numneu;
  size_t state_bytes = (size_t)n_paths * params_.ne * params_.nrhos * su_size * sizeof(double);

  if (d_states_) cudaFree(d_states_);
  NUSQUIDS_CUDA_CHECK(cudaMalloc(&d_states_, state_bytes));

  if (d_paths_) cudaFree(d_paths_);
  NUSQUIDS_CUDA_CHECK(cudaMalloc(&d_paths_, n_paths * sizeof(PathDeviceData)));
}

void Propagator::freeBatch() {
  if (d_states_) { cudaFree(d_states_); d_states_ = nullptr; }
  if (d_paths_) { cudaFree(d_paths_); d_paths_ = nullptr; }
}

void Propagator::uploadProfiles(const std::vector<GPUDensityProfile>& profiles,
                                std::vector<PathDeviceData>& host_paths) {
  host_paths.resize(profiles.size());

  // Calculate total spline data needed
  size_t total_spline_doubles = 0;
  for (const auto& prof : profiles) {
    if (prof.type == ProfileType::TABULATED) {
      // density: x + 4 coefficients = 5 arrays
      // ye: same = 5 arrays (shares x with density)
      // Each array has n nodes
      total_spline_doubles += prof.density_spline.n * 9; // x + 4*density + 4*ye
      for (int t = 0; t < prof.n_targets; t++) {
        if (!prof.target_splines.empty())
          total_spline_doubles += prof.density_spline.n * 4; // 4 coefficients per target
      }
    }
  }

  // Reallocate spline buffer if needed
  if (total_spline_doubles * sizeof(double) > spline_buffer_size_) {
    if (d_spline_buffer_) cudaFree(d_spline_buffer_);
    spline_buffer_size_ = total_spline_doubles * sizeof(double);
    NUSQUIDS_CUDA_CHECK(cudaMalloc(&d_spline_buffer_, spline_buffer_size_));
  }

  // Upload spline data and build PathDeviceData structs
  size_t offset = 0;
  std::vector<double> spline_host(total_spline_doubles);

  for (size_t i = 0; i < profiles.size(); i++) {
    const auto& prof = profiles[i];
    auto& path = host_paths[i];

    path.profile.type = prof.type;
    path.profile.constant_density = prof.constant_density;
    path.profile.constant_ye = prof.constant_ye;
    path.profile.n_targets = prof.n_targets;
    path.xini = prof.xini;
    path.xend = prof.xend;
    path.time_offset = 0.0;

    for (int t = 0; t < MAX_TARGETS; t++) {
      path.profile.constant_target_fractions[t] =
        (t < (int)prof.constant_target_fractions.size()) ?
        prof.constant_target_fractions[t] : 0.0;
      path.profile.target_a0[t] = nullptr;
      path.profile.target_a1[t] = nullptr;
      path.profile.target_a2[t] = nullptr;
      path.profile.target_a3[t] = nullptr;
    }

    if (prof.type == ProfileType::TABULATED) {
      int n = prof.density_spline.n;
      path.profile.n_nodes = n;
      path.profile.xini = prof.xini;
      path.profile.xend = prof.xend;

      // Copy spline data to staging buffer and record device pointers
      auto copyArray = [&](const std::vector<double>& src) -> const double* {
        const double* dev_ptr = d_spline_buffer_ + offset;
        std::copy(src.begin(), src.end(), spline_host.data() + offset);
        offset += src.size();
        return dev_ptr;
      };

      path.profile.density_x = copyArray(prof.density_spline.x);
      path.profile.density_a0 = copyArray(prof.density_spline.a0);
      path.profile.density_a1 = copyArray(prof.density_spline.a1);
      path.profile.density_a2 = copyArray(prof.density_spline.a2);
      path.profile.density_a3 = copyArray(prof.density_spline.a3);

      path.profile.ye_x = path.profile.density_x; // shared x nodes
      path.profile.ye_a0 = copyArray(prof.ye_spline.a0);
      path.profile.ye_a1 = copyArray(prof.ye_spline.a1);
      path.profile.ye_a2 = copyArray(prof.ye_spline.a2);
      path.profile.ye_a3 = copyArray(prof.ye_spline.a3);

      for (int t = 0; t < prof.n_targets && t < (int)prof.target_splines.size(); t++) {
        path.profile.target_a0[t] = copyArray(prof.target_splines[t].a0);
        path.profile.target_a1[t] = copyArray(prof.target_splines[t].a1);
        path.profile.target_a2[t] = copyArray(prof.target_splines[t].a2);
        path.profile.target_a3[t] = copyArray(prof.target_splines[t].a3);
      }
    }
  }

  // Upload all spline data in one transfer
  if (total_spline_doubles > 0) {
    NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(d_spline_buffer_, spline_host.data(),
                                         total_spline_doubles * sizeof(double),
                                         cudaMemcpyHostToDevice, stream_));
  }

  // Upload PathDeviceData array
  NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(d_paths_, host_paths.data(),
                                       host_paths.size() * sizeof(PathDeviceData),
                                       cudaMemcpyHostToDevice, stream_));
}

void Propagator::evolveBatch(const std::vector<GPUDensityProfile>& profiles,
                             const double* initial_states,
                             double* final_states,
                             int n_paths) {
  NUSQUIDS_CUDA_CHECK(cudaSetDevice(device_id_));

  int su_size = params_.numneu * params_.numneu;
  size_t state_bytes = (size_t)n_paths * params_.ne * params_.nrhos * su_size * sizeof(double);

  // Allocate device memory for this batch
  allocateBatch(n_paths);

  // Upload profiles
  std::vector<PathDeviceData> host_paths;
  uploadProfiles(profiles, host_paths);

  // Upload initial states
  NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(d_states_, initial_states, state_bytes,
                                       cudaMemcpyHostToDevice, stream_));

  // Synchronize before kernel launch
  NUSQUIDS_CUDA_CHECK(cudaStreamSynchronize(stream_));

  // Launch evolution kernel
  // Pass interaction data if cross-sections have been uploaded (n_targets > 0)
  const InteractionDataGPU* int_data_ptr =
    (interaction_data_.n_targets > 0) ? &interaction_data_ : nullptr;
  launchEvolve(params_, d_paths_, d_H0_array_, d_b1_proj_,
               int_data_ptr,
               solver_config_, d_states_, n_paths, params_.numneu, stream_);

  // Download final states
  NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(final_states, d_states_, state_bytes,
                                       cudaMemcpyDeviceToHost, stream_));

  // Wait for completion
  NUSQUIDS_CUDA_CHECK(cudaStreamSynchronize(stream_));
}

}} // namespace nusquids::cuda
