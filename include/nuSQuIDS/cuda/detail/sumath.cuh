/******************************************************************************
*   SU(N) algebra operations on GPU for nuSQuIDS CUDA backend.
*   Hardcoded, fully-unrolled operations for SU(3) and SU(4).
*   Adapted from CUDAnuSQuIDS (GPL-3.0).
*
*   All operations work on double* arrays:
*     SU(3): 9 components (dim^2)
*     SU(4): 16 components (dim^2)
******************************************************************************/

#ifndef NUSQUIDS_CUDA_SUMATH_CUH
#define NUSQUIDS_CUDA_SUMATH_CUH

#include <cuda_runtime.h>

namespace nusquids { namespace cuda {

// ============================================================
// SU(N) representation sizes
// ============================================================

__host__ __device__ __forceinline__
constexpr int suSize(int nflv) { return nflv * nflv; }

// ============================================================
// Basic operations: scale, add, copy
// ============================================================

template<int NCOMP>
__device__ __forceinline__
void suScale(const double* __restrict__ a, double s, double* __restrict__ result) {
  #pragma unroll
  for (int i = 0; i < NCOMP; i++)
    result[i] = a[i] * s;
}

template<int NCOMP>
__device__ __forceinline__
void suAdd(const double* __restrict__ a, const double* __restrict__ b,
           double* __restrict__ result) {
  #pragma unroll
  for (int i = 0; i < NCOMP; i++)
    result[i] = a[i] + b[i];
}

template<int NCOMP>
__device__ __forceinline__
void suAddScaled(const double* __restrict__ a, double sa,
                 const double* __restrict__ b, double sb,
                 double* __restrict__ result) {
  #pragma unroll
  for (int i = 0; i < NCOMP; i++)
    result[i] = sa * a[i] + sb * b[i];
}

template<int NCOMP>
__device__ __forceinline__
void suCopy(const double* __restrict__ src, double* __restrict__ dst) {
  #pragma unroll
  for (int i = 0; i < NCOMP; i++)
    dst[i] = src[i];
}

template<int NCOMP>
__device__ __forceinline__
void suZero(double* __restrict__ a) {
  #pragma unroll
  for (int i = 0; i < NCOMP; i++)
    a[i] = 0.0;
}

// ============================================================
// SU(3) operations (9 components)
// ============================================================

/// Trace Tr(A·B) for SU(3) represented as 9-component vectors.
/// Matches SQuIDS SUTrace: dim*a[0]*b[0] + 2*sum_{i>0} a[i]*b[i]
__device__ __forceinline__
double suTrace3(const double* __restrict__ a, const double* __restrict__ b) {
  double result = 3.0 * a[0] * b[0];
  #pragma unroll
  for (int i = 1; i < 9; i++)
    result += 2.0 * a[i] * b[i];
  return result;
}

/// Trace Tr(A·B) for SU(4) represented as 16-component vectors.
/// Matches SQuIDS SUTrace: dim*a[0]*b[0] + 2*sum_{i>0} a[i]*b[i]
__device__ __forceinline__
double suTrace4(const double* __restrict__ a, const double* __restrict__ b) {
  double result = 4.0 * a[0] * b[0];
  #pragma unroll
  for (int i = 1; i < 16; i++)
    result += 2.0 * a[i] * b[i];
  return result;
}

/// Generic trace dispatcher
template<int NFLV>
__device__ __forceinline__
double suTrace(const double* __restrict__ a, const double* __restrict__ b);

template<>
__device__ __forceinline__
double suTrace<3>(const double* __restrict__ a, const double* __restrict__ b) {
  return suTrace3(a, b);
}

template<>
__device__ __forceinline__
double suTrace<4>(const double* __restrict__ a, const double* __restrict__ b) {
  return suTrace4(a, b);
}

// ============================================================
// Commutator: result = i[A, B]
// Anticommutator: result = {A, B}
// These are the core operations used in the evolution equation.
// ============================================================

/// i[A, B] for SU(3) — commutator (9 components)
/// Structure constants are hardcoded following the SQuIDS convention.
__device__
void suCommutator3(const double* __restrict__ a,
                   const double* __restrict__ b,
                   double* __restrict__ result);

/// {A, B} for SU(3) — anticommutator (9 components)
__device__
void suAnticommutator3(const double* __restrict__ a,
                       const double* __restrict__ b,
                       double* __restrict__ result);

/// i[A, B] for SU(4) — commutator (16 components)
__device__
void suCommutator4(const double* __restrict__ a,
                   const double* __restrict__ b,
                   double* __restrict__ result);

/// {A, B} for SU(4) — anticommutator (16 components)
__device__
void suAnticommutator4(const double* __restrict__ a,
                       const double* __restrict__ b,
                       double* __restrict__ result);

// ============================================================
// Unitary evolution: evolve state with Hamiltonian
// result = exp(-i H t) * state * exp(i H t)
//
// In the SU(N) basis this corresponds to rotating each
// off-diagonal component by the appropriate phase differences.
// ============================================================

/// Evolve an SU(3) vector under a diagonal Hamiltonian
/// \param state  Input SU(3) vector (9 components)
/// \param phases  Pre-computed phase differences from H0 eigenvalues and time
/// \param result Output evolved SU(3) vector
__device__
void suEvolve3(const double* __restrict__ state,
               const double* __restrict__ phases,
               double* __restrict__ result);

/// Evolve an SU(4) vector under a diagonal Hamiltonian
__device__
void suEvolve4(const double* __restrict__ state,
               const double* __restrict__ phases,
               double* __restrict__ result);

}} // namespace nusquids::cuda

#endif // NUSQUIDS_CUDA_SUMATH_CUH
