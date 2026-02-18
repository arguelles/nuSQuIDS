/******************************************************************************
*   GPU body profiles for nuSQuIDS CUDA backend.
*   Pre-tabulates density/ye/composition profiles on CPU using Akima splines,
*   then evaluates on GPU. This approach supports ANY body type without
*   device-side virtual dispatch.
*   Adapted from CUDAnuSQuIDS EarthAtm implementation (GPL-3.0).
******************************************************************************/

#ifndef NUSQUIDS_CUDA_BODY_GPU_CUH
#define NUSQUIDS_CUDA_BODY_GPU_CUH

#include <cuda_runtime.h>
#include <vector>
#include <map>

namespace nusquids { namespace cuda {

// ============================================================
// CPU-side: profile extraction and Akima spline fitting
// ============================================================

/// Type of density profile detected
enum class ProfileType {
  VACUUM,     ///< All density samples are zero
  CONSTANT,   ///< Density and ye are constant along the track
  TABULATED   ///< General case — use Akima spline interpolation
};

/// Akima spline coefficients for one quantity along a track
struct AkimaCoeffs {
  std::vector<double> x;    ///< Node positions
  std::vector<double> a0;   ///< Coefficient a0 (value at node)
  std::vector<double> a1;   ///< Coefficient a1 (first derivative)
  std::vector<double> a2;   ///< Coefficient a2 (second derivative / 2)
  std::vector<double> a3;   ///< Coefficient a3 (third derivative / 6)
  int n;                     ///< Number of nodes
};

/// CPU-side representation of a density profile along one track
struct GPUDensityProfile {
  ProfileType type;
  double constant_density;   ///< Used for CONSTANT type
  double constant_ye;        ///< Used for CONSTANT type

  AkimaCoeffs density_spline; ///< Akima spline for density
  AkimaCoeffs ye_spline;      ///< Akima spline for ye

  /// Per-target composition fractions (keyed by target index)
  /// Each target gets its own Akima spline if composition varies
  std::vector<double> constant_target_fractions;
  std::vector<AkimaCoeffs> target_splines;
  int n_targets;

  double xini;  ///< Track start position
  double xend;  ///< Track end position
};

/// Fit Akima spline to (x, y) data
AkimaCoeffs fitAkimaSpline(const std::vector<double>& x,
                            const std::vector<double>& y);

/// Extract density profile from a Body+Track by sampling
/// \param body The body to sample
/// \param track The track to sample along
/// \param n_samples Number of sample points (default 1000)
/// \param target_indices Indices into int_struct.targets for composition sampling
GPUDensityProfile extractProfile(const void* body_ptr, const void* track_ptr,
                                 int n_samples = 1000);

// ============================================================
// GPU-side: compact profile data and evaluation
// ============================================================

/// Maximum number of spline nodes (set generously; most profiles need ~100-1000)
static constexpr int MAX_PROFILE_NODES = 2048;
/// Maximum number of nuclear targets
static constexpr int MAX_TARGETS = 16;

/// Device-side compact profile for one track
struct GPUDensityProfileDevice {
  ProfileType type;

  // For CONSTANT type
  double constant_density;
  double constant_ye;

  // For TABULATED type: Akima spline data (stored contiguously)
  int n_nodes;
  double xini, xend;

  // Pointers into device memory (set during upload)
  const double* density_x;
  const double* density_a0;
  const double* density_a1;
  const double* density_a2;
  const double* density_a3;

  const double* ye_x;
  const double* ye_a0;
  const double* ye_a1;
  const double* ye_a2;
  const double* ye_a3;

  // Target composition
  int n_targets;
  double constant_target_fractions[MAX_TARGETS];
  const double* target_a0[MAX_TARGETS];
  const double* target_a1[MAX_TARGETS];
  const double* target_a2[MAX_TARGETS];
  const double* target_a3[MAX_TARGETS];
};

/// Evaluate Akima spline at position x
__device__ __forceinline__
double evaluateAkimaSpline(double x, int n_nodes,
                           const double* __restrict__ xs,
                           const double* __restrict__ a0,
                           const double* __restrict__ a1,
                           const double* __restrict__ a2,
                           const double* __restrict__ a3) {
  // Binary search for interval
  int lo = 0, hi = n_nodes - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (xs[mid] <= x)
      lo = mid;
    else
      hi = mid;
  }
  double dx = x - xs[lo];
  return a0[lo] + dx * (a1[lo] + dx * (a2[lo] + dx * a3[lo]));
}

/// Evaluate density at position x
__device__ __forceinline__
double evaluateDensity(const GPUDensityProfileDevice& profile, double x) {
  if (profile.type == ProfileType::VACUUM)
    return 0.0;
  if (profile.type == ProfileType::CONSTANT)
    return profile.constant_density;
  return evaluateAkimaSpline(x, profile.n_nodes,
                              profile.density_x, profile.density_a0,
                              profile.density_a1, profile.density_a2,
                              profile.density_a3);
}

/// Evaluate electron fraction (Ye) at position x
__device__ __forceinline__
double evaluateYe(const GPUDensityProfileDevice& profile, double x) {
  if (profile.type == ProfileType::VACUUM)
    return 0.0;
  if (profile.type == ProfileType::CONSTANT)
    return profile.constant_ye;
  return evaluateAkimaSpline(x, profile.n_nodes,
                              profile.ye_x, profile.ye_a0,
                              profile.ye_a1, profile.ye_a2,
                              profile.ye_a3);
}

/// Evaluate target fraction at position x
__device__ __forceinline__
double evaluateTargetFraction(const GPUDensityProfileDevice& profile,
                              int target, double x) {
  if (profile.type != ProfileType::TABULATED || target >= profile.n_targets)
    return profile.constant_target_fractions[target];
  if (profile.target_a0[target] == nullptr)
    return profile.constant_target_fractions[target];
  return evaluateAkimaSpline(x, profile.n_nodes,
                              profile.density_x, // reuse x nodes
                              profile.target_a0[target],
                              profile.target_a1[target],
                              profile.target_a2[target],
                              profile.target_a3[target]);
}

}} // namespace nusquids::cuda

#endif // NUSQUIDS_CUDA_BODY_GPU_CUH
