/******************************************************************************
*   GPU memory management utilities for nuSQuIDS CUDA backend.
*   RAII wrappers for device and pinned host memory.
*   Adapted from CUDAnuSQuIDS (GPL-3.0).
******************************************************************************/

#ifndef NUSQUIDS_CUDA_MEMORY_CUH
#define NUSQUIDS_CUDA_MEMORY_CUH

#include <cuda_runtime.h>
#include <memory>
#include <stdexcept>
#include <string>

namespace nusquids { namespace cuda {

/// Check CUDA error and throw on failure
inline void checkCuda(cudaError_t err, const char* file, int line) {
  if (err != cudaSuccess) {
    throw std::runtime_error(
      std::string("CUDA error at ") + file + ":" + std::to_string(line) +
      ": " + cudaGetErrorString(err));
  }
}

#define NUSQUIDS_CUDA_CHECK(err) nusquids::cuda::checkCuda(err, __FILE__, __LINE__)

/// Custom deleter for cudaFree
struct DeviceDeleter {
  void operator()(void* p) const {
    if (p) cudaFree(p);
  }
};

/// Custom deleter for cudaFreeHost
struct PinnedDeleter {
  void operator()(void* p) const {
    if (p) cudaFreeHost(p);
  }
};

/// RAII wrapper for device memory
template<typename T>
using unique_dev_ptr = std::unique_ptr<T[], DeviceDeleter>;

/// RAII wrapper for pinned host memory
template<typename T>
using unique_pinned_ptr = std::unique_ptr<T[], PinnedDeleter>;

/// Allocate device memory
template<typename T>
unique_dev_ptr<T> make_device(size_t count) {
  T* ptr = nullptr;
  NUSQUIDS_CUDA_CHECK(cudaMalloc(&ptr, count * sizeof(T)));
  return unique_dev_ptr<T>(ptr);
}

/// Allocate pinned (page-locked) host memory
template<typename T>
unique_pinned_ptr<T> make_pinned(size_t count) {
  T* ptr = nullptr;
  NUSQUIDS_CUDA_CHECK(cudaMallocHost(&ptr, count * sizeof(T)));
  return unique_pinned_ptr<T>(ptr);
}

/// Pitched 2D allocation info
struct PitchedAlloc {
  void* ptr;
  size_t pitch;   ///< Pitch in bytes
  size_t width;   ///< Width in bytes
  size_t height;
};

/// Allocate pitched 2D device memory for coalesced access
inline PitchedAlloc make_device_pitched(size_t width_bytes, size_t height) {
  PitchedAlloc result;
  result.width = width_bytes;
  result.height = height;
  NUSQUIDS_CUDA_CHECK(cudaMallocPitch(&result.ptr, &result.pitch,
                                       width_bytes, height));
  return result;
}

/// Free pitched allocation
inline void free_pitched(PitchedAlloc& alloc) {
  if (alloc.ptr) {
    cudaFree(alloc.ptr);
    alloc.ptr = nullptr;
  }
}

/// Copy host→device
template<typename T>
void copyToDevice(T* dst, const T* src, size_t count, cudaStream_t stream = 0) {
  NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(dst, src, count * sizeof(T),
                                        cudaMemcpyHostToDevice, stream));
}

/// Copy device→host
template<typename T>
void copyToHost(T* dst, const T* src, size_t count, cudaStream_t stream = 0) {
  NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(dst, src, count * sizeof(T),
                                        cudaMemcpyDeviceToHost, stream));
}

}} // namespace nusquids::cuda

#endif // NUSQUIDS_CUDA_MEMORY_CUH
