/******************************************************************************
*   CUDABackend implementation for nuSQuIDS.
*   GPU detection, multi-GPU dispatch, data marshaling.
*   State extraction and writeback are handled by the caller (nuSQUIDSAtm).
******************************************************************************/

#include "nuSQuIDS/cuda/cuda_backend.h"
#include "nuSQuIDS/cuda/detail/propagator.cuh"
#include "nuSQuIDS/cuda/detail/body_gpu.cuh"
#include "nuSQuIDS/cuda/detail/memory.cuh"

#include <cuda_runtime.h>
#include <thread>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <sstream>

namespace nusquids {

// ============================================================
// Akima spline fitting (CPU-side)
// ============================================================

namespace cuda {

AkimaCoeffs fitAkimaSpline(const std::vector<double>& x,
                            const std::vector<double>& y) {
  int n = (int)x.size();
  AkimaCoeffs result;
  result.n = n;
  result.x = x;
  result.a0 = y;
  result.a1.resize(n, 0.0);
  result.a2.resize(n, 0.0);
  result.a3.resize(n, 0.0);

  if (n < 2) return result;
  if (n == 2) {
    result.a1[0] = (y[1] - y[0]) / (x[1] - x[0]);
    return result;
  }

  // Compute slopes between consecutive points
  std::vector<double> m(n + 3);
  for (int i = 0; i < n - 1; i++)
    m[i + 2] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);

  // Extrapolate end slopes
  m[1] = 2.0 * m[2] - m[3];
  m[0] = 2.0 * m[1] - m[2];
  m[n + 1] = 2.0 * m[n] - m[n - 1];
  m[n + 2] = 2.0 * m[n + 1] - m[n];

  // Compute Akima weights and derivatives
  for (int i = 0; i < n; i++) {
    double w1 = fabs(m[i + 3] - m[i + 2]);
    double w2 = fabs(m[i + 1] - m[i]);
    if (w1 + w2 > 1e-30)
      result.a1[i] = (w1 * m[i + 1] + w2 * m[i + 2]) / (w1 + w2);
    else
      result.a1[i] = 0.5 * (m[i + 1] + m[i + 2]);
  }

  // Compute cubic coefficients
  for (int i = 0; i < n - 1; i++) {
    double dx = x[i + 1] - x[i];
    double dy = y[i + 1] - y[i];
    result.a2[i] = (3.0 * dy / dx - 2.0 * result.a1[i] - result.a1[i + 1]) / dx;
    result.a3[i] = (result.a1[i] + result.a1[i + 1] - 2.0 * dy / dx) / (dx * dx);
  }

  return result;
}

} // namespace cuda

// ============================================================
// CUDABackend::Impl
// ============================================================

struct CUDABackend::Impl {
  Config config;
  std::vector<int> device_ids;
  std::vector<std::unique_ptr<cuda::Propagator>> propagators;

  Impl(Config config) : config(config) {
    // Detect available GPUs
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0)
      throw std::runtime_error("CUDABackend: No CUDA-capable GPU detected");

    if (config.device_ids.empty()) {
      device_ids.resize(device_count);
      std::iota(device_ids.begin(), device_ids.end(), 0);
    } else {
      device_ids = config.device_ids;
      for (int id : device_ids) {
        if (id >= device_count)
          throw std::runtime_error("CUDABackend: Invalid device ID " + std::to_string(id));
      }
    }

    for (int id : device_ids) {
      propagators.push_back(
        std::make_unique<cuda::Propagator>(id, config.batch_size_limit));
    }
  }

  /// Convert GPUPathData to cuda::GPUDensityProfile for the propagator
  cuda::GPUDensityProfile convertPath(const GPUPathData& path) {
    cuda::GPUDensityProfile profile;
    profile.xini = path.xini;
    profile.xend = path.xend;

    if (path.xend <= path.xini || path.n_density_samples == 0) {
      profile.type = cuda::ProfileType::VACUUM;
      profile.constant_density = 0.0;
      profile.constant_ye = 0.0;
      profile.n_targets = 0;
      return profile;
    }

    // Check for vacuum (all densities zero)
    bool all_zero = true;
    for (int i = 0; i < path.n_density_samples; i++) {
      if (path.density_vals[i] > 1e-30) { all_zero = false; break; }
    }
    if (all_zero) {
      profile.type = cuda::ProfileType::VACUUM;
      profile.constant_density = 0.0;
      profile.constant_ye = 0.0;
      profile.n_targets = 0;
      return profile;
    }

    // Check for constant density
    bool all_same_density = true;
    bool all_same_ye = true;
    for (int i = 1; i < path.n_density_samples; i++) {
      if (fabs(path.density_vals[i] - path.density_vals[0]) > 1e-10 * fabs(path.density_vals[0]))
        all_same_density = false;
      if (fabs(path.ye_vals[i] - path.ye_vals[0]) > 1e-10)
        all_same_ye = false;
    }
    if (all_same_density && all_same_ye) {
      profile.type = cuda::ProfileType::CONSTANT;
      profile.constant_density = path.density_vals[0];
      profile.constant_ye = path.ye_vals[0];
      profile.n_targets = 0;
      return profile;
    }

    // General case: fit Akima splines
    profile.type = cuda::ProfileType::TABULATED;
    profile.constant_density = 0.0;
    profile.constant_ye = 0.0;
    profile.density_spline = cuda::fitAkimaSpline(path.density_x, path.density_vals);
    profile.ye_spline = cuda::fitAkimaSpline(path.density_x, path.ye_vals);
    profile.n_targets = 0;

    return profile;
  }
};

// ============================================================
// CUDABackend public API
// ============================================================

CUDABackend::CUDABackend(Config config)
  : impl_(std::make_unique<Impl>(config)) {}

CUDABackend::~CUDABackend() = default;

CUDABackend::CUDABackend(CUDABackend&&) noexcept = default;
CUDABackend& CUDABackend::operator=(CUDABackend&&) noexcept = default;

void CUDABackend::Evolve(double* states,
                          const std::vector<GPUPathData>& paths,
                          const double* H0_array,
                          const double* b1_proj,
                          int n_paths, int ne, int nrhos, int numneu,
                          double HI_constants, int NT_type,
                          double rel_error, double abs_error) {
  if (n_paths == 0) return;

  int su_size = numneu * numneu;
  int state_per_path = nrhos * ne * su_size;

  // Build PhysicsParams
  cuda::PhysicsParams params;
  params.numneu = numneu;
  params.nrhos = nrhos;
  params.ne = ne;
  params.su_size = su_size;
  params.HI_constants = HI_constants;
  params.NT_type = NT_type;
  params.iinteraction = false;
  params.ioscillations = true;
  params.tauregeneration = false;
  params.iglashow = false;
  params.basis = 2; // interaction basis

  // Build solver config — use caller-provided tolerances if given
  cuda::SolverConfig solver_config;
  solver_config.h_initial = 1.0;
  solver_config.h_min = 1e-10;
  solver_config.h_max = 1e20;
  solver_config.rel_error = (rel_error > 0.0) ? rel_error : impl_->config.rel_error;
  solver_config.abs_error = (abs_error > 0.0) ? abs_error : impl_->config.abs_error;
  solver_config.max_steps = 1000000;

  // Convert paths to GPU density profiles
  std::vector<cuda::GPUDensityProfile> profiles(n_paths);
  for (int i = 0; i < n_paths; i++)
    profiles[i] = impl_->convertPath(paths[i]);

  // Upload shared data to all GPUs
  cuda::InteractionDataGPU empty_interaction;
  memset(&empty_interaction, 0, sizeof(empty_interaction));
  for (auto& prop : impl_->propagators) {
    prop->uploadSharedData(params, empty_interaction,
                          H0_array, b1_proj,
                          solver_config);
  }

  int n_gpus = impl_->propagators.size();

  // Allocate output buffer
  std::vector<double> final_states(n_paths * state_per_path);

  if (n_gpus == 1) {
    auto& prop = impl_->propagators[0];
    int batch_limit = impl_->config.batch_size_limit;

    for (int batch_start = 0; batch_start < n_paths; batch_start += batch_limit) {
      int batch_size = std::min(batch_limit, n_paths - batch_start);
      std::vector<cuda::GPUDensityProfile> batch_profiles(
        profiles.begin() + batch_start,
        profiles.begin() + batch_start + batch_size);

      prop->evolveBatch(batch_profiles,
                        states + batch_start * state_per_path,
                        final_states.data() + batch_start * state_per_path,
                        batch_size);
    }
  } else {
    // Multi-GPU: distribute via threads
    std::vector<std::thread> gpu_threads;

    for (int gpu = 0; gpu < n_gpus; gpu++) {
      gpu_threads.emplace_back([&, gpu]() {
        auto& prop = impl_->propagators[gpu];

        // Round-robin path assignment
        std::vector<int> my_indices;
        for (int i = gpu; i < n_paths; i += n_gpus)
          my_indices.push_back(i);

        if (my_indices.empty()) return;

        int my_n = my_indices.size();
        std::vector<cuda::GPUDensityProfile> my_profiles(my_n);
        std::vector<double> my_states(my_n * state_per_path);

        for (int j = 0; j < my_n; j++) {
          my_profiles[j] = profiles[my_indices[j]];
          std::copy(states + my_indices[j] * state_per_path,
                    states + (my_indices[j] + 1) * state_per_path,
                    my_states.begin() + j * state_per_path);
        }

        std::vector<double> my_final(my_n * state_per_path);
        prop->evolveBatch(my_profiles, my_states.data(),
                         my_final.data(), my_n);

        for (int j = 0; j < my_n; j++) {
          std::copy(my_final.begin() + j * state_per_path,
                    my_final.begin() + (j + 1) * state_per_path,
                    final_states.begin() + my_indices[j] * state_per_path);
        }
      });
    }

    for (auto& t : gpu_threads)
      t.join();
  }

  // Copy evolved states back to the caller's buffer
  std::copy(final_states.begin(), final_states.end(), states);
}

bool CUDABackend::IsAvailable() {
  int count = 0;
  cudaError_t err = cudaGetDeviceCount(&count);
  return (err == cudaSuccess && count > 0);
}

int CUDABackend::DeviceCount() {
  int count = 0;
  cudaGetDeviceCount(&count);
  return count;
}

std::string CUDABackend::DeviceInfo() {
  int count = 0;
  cudaGetDeviceCount(&count);
  std::ostringstream oss;
  for (int i = 0; i < count; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    oss << "GPU " << i << ": " << prop.name
        << " (SM " << prop.major << "." << prop.minor
        << ", " << (prop.totalGlobalMem >> 20) << " MB)\n";
  }
  return oss.str();
}

} // namespace nusquids
