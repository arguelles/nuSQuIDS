/******************************************************************************
*   CUDA kernel implementations for nuSQuIDS.
*   Interaction picture evolution with correct SU(3) algebra.
*   Formulas directly from SQuIDS source (iConmutatorSU3.txt, SinCosEvolSU3.txt).
*   Adapted from CUDAnuSQuIDS (GPL-3.0).
******************************************************************************/

#include "nuSQuIDS/cuda/detail/kernels.cuh"
#include "nuSQuIDS/cuda/detail/sumath.cuh"
#include "nuSQuIDS/cuda/detail/solver_gpu.cuh"
#include "nuSQuIDS/cuda/detail/physics.cuh"
#include "nuSQuIDS/cuda/detail/memory.cuh"
#include "nuSQuIDS/cuda/detail/interactions_gpu.cuh"

namespace nusquids { namespace cuda {

// ============================================================
// Constants
// ============================================================

// ============================================================
// iCommutator for SU(3) — exact from SQuIDS iConmutatorSU3.txt
// Computes result = i[A, B] where A=state, B=Hamiltonian
// Convention matches SQuIDS: dstate = iCommutator(estate, HI)
// ============================================================

__device__ __forceinline__
void iCommutatorSU3(const double* __restrict__ A,
                    const double* __restrict__ B,
                    double* __restrict__ R)
{
  const double s3 = 1.7320508075688772;

  R[0] = 0.0;

  R[1] = A[7]*B[2] + 2.0*A[4]*B[3] - 2.0*A[3]*B[4]
       + A[6]*B[5] - A[5]*B[6] - A[2]*B[7];

  R[2] = -A[7]*B[1] - A[5]*B[3] + A[3]*B[5]
       + (A[4] + s3*A[8])*B[6] + A[1]*B[7]
       - A[6]*(B[4] + s3*B[8]);

  R[3] = -2.0*A[4]*B[1] + A[5]*B[2] + 2.0*A[1]*B[4]
       - A[2]*B[5] + A[7]*B[6] - A[6]*B[7];

  R[4] = 2.0*A[3]*B[1] + A[6]*B[2] - 2.0*A[1]*B[3]
       - A[7]*B[5] - A[2]*B[6] + A[5]*B[7];

  R[5] = -A[6]*B[1] - A[3]*B[2] + A[2]*B[3]
       + A[1]*B[6] + (-A[4] + s3*A[8])*B[7]
       + A[7]*(B[4] - s3*B[8]);

  R[6] = A[5]*B[1] - (A[4] + s3*A[8])*B[2]
       - A[7]*B[3] - A[1]*B[5] + A[3]*B[7]
       + A[2]*(B[4] + s3*B[8]);

  R[7] = A[2]*B[1] - A[1]*B[2] + A[6]*B[3]
       + (A[4] - s3*A[8])*B[5] - A[3]*B[6]
       + A[5]*(-B[4] + s3*B[8]);

  R[8] = s3*(A[6]*B[2] + A[7]*B[5]
       - A[2]*B[6] - A[5]*B[7]);
}

// ============================================================
// Anticommutator {A, B} for SU(3) — from SQuIDS AnticonmutatorSU3.txt
// Used for absorption: dρ/dt -= {Γ, ρ}
// ============================================================

__device__ __forceinline__
void antiCommutatorSU3(const double* __restrict__ A,
                       const double* __restrict__ B,
                       double* __restrict__ R)
{
  const double s3 = 1.7320508075688772;

  R[0] = (2.0*(3.0*A[0]*B[0] + 2.0*(A[1]*B[1] + A[2]*B[2] + A[3]*B[3]
        + A[4]*B[4] + A[5]*B[5] + A[6]*B[6] + A[7]*B[7] + A[8]*B[8])))/3.0;

  R[1] = 2.0*A[0]*B[1] + (2.0*A[8]*B[1])/s3 + A[5]*B[2] + A[2]*B[5]
       + A[7]*B[6] + A[6]*B[7] + A[1]*(2.0*B[0] + (2.0*B[8])/s3);

  R[2] = A[5]*B[1] + ((6.0*A[0] + 3.0*A[4] - s3*A[8])*B[2])/3.0
       - A[7]*B[3] + A[1]*B[5] - A[3]*B[7]
       + A[2]*(2.0*B[0] + B[4] - B[8]/s3);

  R[3] = -(A[7]*B[2]) + 2.0*A[0]*B[3] + (2.0*A[8]*B[3])/s3
       + A[6]*B[5] + A[5]*B[6] - A[2]*B[7]
       + A[3]*(2.0*B[0] + (2.0*B[8])/s3);

  R[4] = A[2]*B[2] + 2.0*A[0]*B[4] + (2.0*A[8]*B[4])/s3
       - A[5]*B[5] + A[6]*B[6] - A[7]*B[7]
       + A[4]*(2.0*B[0] + (2.0*B[8])/s3);

  R[5] = A[2]*B[1] + A[1]*B[2] + A[6]*B[3]
       + ((6.0*A[0] - 3.0*A[4] - s3*A[8])*B[5])/3.0 + A[3]*B[6]
       + A[5]*(2.0*B[0] - B[4] - B[8]/s3);

  R[6] = A[7]*B[1] + A[5]*B[3] + A[3]*B[5]
       + ((6.0*A[0] + 3.0*A[4] - s3*A[8])*B[6])/3.0 + A[1]*B[7]
       + A[6]*(2.0*B[0] + B[4] - B[8]/s3);

  R[7] = A[6]*B[1] - A[3]*B[2] - A[2]*B[3] + A[1]*B[6]
       + ((6.0*A[0] - 3.0*A[4] - s3*A[8])*B[7])/3.0
       + A[7]*(2.0*B[0] - B[4] - B[8]/s3);

  R[8] = (2.0*A[1]*B[1] - A[2]*B[2] + 2.0*A[3]*B[3] + 2.0*A[4]*B[4]
       - A[5]*B[5] - A[6]*B[6] - A[7]*B[7]
       + 2.0*A[8]*(s3*B[0] - B[8]) + 2.0*s3*A[0]*B[8]) / s3;
}

// ============================================================
// Compute evolved projectors at position x_eval for SU(3)
// Returns the SinCosEvol-rotated projectors in evol_proj_out
// ============================================================

__device__ __forceinline__
void evolveProjectorsSU3(double x_eval, double xini,
                         const double* __restrict__ H0,
                         const double* __restrict__ b1_proj,
                         int numneu,
                         double* __restrict__ evol_proj_out)
{
  const double s3 = 1.7320508075688772;
  double dt = x_eval - xini;
  double w12 = 2.0 * H0[4];
  double w13 = H0[4] + s3 * H0[8];
  double w23 = H0[4] - s3 * H0[8];

  double CX0, SX0, CX1, SX1, CX2, SX2;
  sincos(w12 * dt, &SX0, &CX0);
  sincos(w13 * dt, &SX1, &CX1);
  sincos(w23 * dt, &SX2, &CX2);

  for (int flv = 0; flv < numneu; flv++) {
    const double* p = b1_proj + flv * 9;
    double* e = evol_proj_out + flv * 9;
    e[0] = p[0];
    e[1] = CX0*p[1] + SX0*p[3];
    e[2] = CX1*p[2] + SX1*p[6];
    e[3] = CX0*p[3] - SX0*p[1];
    e[4] = p[4];
    e[5] = CX2*p[5] - SX2*p[7];
    e[6] = CX1*p[6] - SX1*p[2];
    e[7] = CX2*p[7] + SX2*p[5];
    e[8] = p[8];
  }
}

// ============================================================
// Compute interaction-picture HI at position x_eval for SU(3)
//
// Evolves flavor projectors using SinCosEvolSU3 formulas, then
// constructs HI = (CC+NC)*evol_proj_e + NC*evol_proj_mu + NC*evol_proj_tau
//
// From nuSQUIDS::HI():
//   CC = HI_constants * density * ye
//   NC = CC * (-0.5 * (1 - ye) / ye)
//   For antineutrinos: CC, NC sign-flipped
// ============================================================

__device__ __forceinline__
void computeHI_SU3(double x_eval, double xini,
                   const double* __restrict__ H0,
                   const double* __restrict__ b1_proj,
                   const GPUDensityProfileDevice& profile,
                   double HI_constants, bool is_antinu,
                   double* __restrict__ HI)
{
  double density = evaluateDensity(profile, x_eval);
  double ye = evaluateYe(profile, x_eval);

  if (density <= 1.0e-30) {
    #pragma unroll
    for (int c = 0; c < 9; c++) HI[c] = 0.0;
    return;
  }

  double CC = HI_constants * density * ye;
  double NC;
  if (ye < 1.0e-10)
    NC = HI_constants * density;
  else
    NC = CC * (-0.5 * (1.0 - ye) / ye);

  if (is_antinu) { CC = -CC; NC = -NC; }

  double weights[3] = {CC + NC, NC, NC};

  // Compute sin/cos for projector evolution (PreSinCosEvolSU3)
  double dt = x_eval - xini;
  const double s3 = 1.7320508075688772;
  double w12 = 2.0 * H0[4];
  double w13 = H0[4] + s3 * H0[8];
  double w23 = H0[4] - s3 * H0[8];

  double CX0, SX0, CX1, SX1, CX2, SX2;
  sincos(w12 * dt, &SX0, &CX0);
  sincos(w13 * dt, &SX1, &CX1);
  sincos(w23 * dt, &SX2, &CX2);

  // Initialize HI = 0
  #pragma unroll
  for (int c = 0; c < 9; c++) HI[c] = 0.0;

  // Accumulate HI = sum_flv weight[flv] * SinCosEvolSU3(b1_proj[flv])
  double evol[9];
  #pragma unroll
  for (int flv = 0; flv < 3; flv++) {
    const double* p = b1_proj + flv * 9;
    double w = weights[flv];

    // SinCosEvolSU3 from SQuIDS
    evol[0] = p[0];
    evol[1] = CX0*p[1] + SX0*p[3];
    evol[2] = CX1*p[2] + SX1*p[6];
    evol[3] = CX0*p[3] - SX0*p[1];
    evol[4] = p[4];
    evol[5] = CX2*p[5] - SX2*p[7];
    evol[6] = CX1*p[6] - SX1*p[2];
    evol[7] = CX2*p[7] + SX2*p[5];
    evol[8] = p[8];

    #pragma unroll
    for (int c = 0; c < 9; c++)
      HI[c] += w * evol[c];
  }
}

// ============================================================
// Electron number density for Glashow (natural units, eV^3).
//
// All branches are precomputed on the CPU side (see the GPU path builder
// in nuSQuIDS.cpp, which mirrors UpdateInteractions and handles isoscalar,
// p/n, body-composition, and nuclear-XS cases in squids::Const natural
// units) and uploaded as either a constant or an Akima spline along the
// track. CPU is the oracle.
//
// Falls back to the legacy isoscalar / p-n derivation from target number
// densities when num_e was not provided (e.g. older callers or test paths
// that skip the precomputation).
// ============================================================

__device__ __forceinline__
double electronNumberDensitySU3(const GPUDensityProfileDevice& profile,
                                double x_eval,
                                const double* __restrict__ target_ndens,
                                int n_targets, double ye)
{
  if (profile.has_num_e)
    return evaluateNumE(profile, x_eval);
  if (n_targets <= 0) return 0.0;
  if (n_targets == 1) return target_ndens[0] * ye;
  return target_ndens[0]; // proton density; in neutral matter n_e = n_p
}

// ============================================================
// Compute inverse interaction lengths from cross sections × density
// invlen[flv] = sum_target density * Na * sigma[target][flv] * target_fraction[target]
//
// For nuSQuIDS, invlen_INT is stored as invlen_CC + invlen_NC.
// HI_constants already encodes sqrt(2)*GF*Na/cm^3 in eV, so we need
// the density in g/cm^3 and sigma in eV^-2. The product is in eV.
//
// When Glashow is enabled and this is the electron antineutrino rho,
// invlen_GR = sigma_GR(e) * n_e is added into invlen[0] (electron), matching
// the CPU in src/nuSQuIDS.cpp:794-803.
// ============================================================

__device__ __forceinline__
void computeInvlenSU3(int ie, int rho, int ne,
                      const double* __restrict__ target_ndens, // [n_targets] precomputed number densities
                      const InteractionDataGPU& idata,
                      bool iglashow, int NT_type, double ye,
                      double num_e,                            // precomputed electron number density (eV^3)
                      double* __restrict__ invlen_out)  // [numneu] output: invlen_INT per flavor
{
  // invlen_INT[flv] = sum over targets: ndens[t] * (sigma_CC + sigma_NC)
  // target_ndens are precomputed on CPU in natural units (eV³) using squids::Const
  for (int flv = 0; flv < 3; flv++) {
    double invlen = 0.0;
    for (int t = 0; t < idata.n_targets; t++) {
      size_t idx = sigma_index(t, rho, flv, ie, idata.nrhos, 3, idata.rounded_ne);
      double sig = idata.d_sigma_CC[idx] + idata.d_sigma_NC[idx];
      invlen += target_ndens[t] * sig;
    }
    invlen_out[flv] = invlen;
  }

  // Glashow adds to electron antineutrino absorption only.
  // NT_type==3 means "both" (rho=0 nu, rho=1 nubar); NT_type==2 means pure antineutrino.
  if (iglashow && idata.d_sigma_GR != nullptr) {
    int gr_rho = (NT_type == 3) ? 1 : 0;
    if (rho == gr_rho && (NT_type == 2 || NT_type == 3)) {
      invlen_out[0] += idata.d_sigma_GR[ie] * num_e;
    }
  }
}

// ============================================================
// Compute GammaRho (absorption term) for SU(3)
// GammaRho = sum_flv 0.5 * invlen_INT[flv] * evol_proj[flv]
//
// The factor 0.5 comes from the master equation: dρ/dt -= {Γ, ρ}
// where the anticommutator doubles the effect.
// ============================================================

__device__ __forceinline__
void computeGammaRhoSU3(const double* __restrict__ invlen,  // [3] invlen_INT per flavor
                         const double* __restrict__ evol_proj,  // [3][9] evolved projectors
                         double* __restrict__ Gamma)  // [9] output
{
  #pragma unroll
  for (int c = 0; c < 9; c++) Gamma[c] = 0.0;

  for (int flv = 0; flv < 3; flv++) {
    double w = 0.5 * invlen[flv];
    const double* ep = evol_proj + flv * 9;
    #pragma unroll
    for (int c = 0; c < 9; c++)
      Gamma[c] += w * ep[c];
  }
}

// ============================================================
// RK4 step for SU(3) interaction picture evolution
// d/dt rho_tilde = iCommutator(rho_tilde, HI_tilde(t))
//
// Optimization: HI at x+h/2 is computed once and reused for k2,k3
// (3 HI evaluations per step instead of 4)
// ============================================================

__device__
void rk4StepSU3(const double* __restrict__ y, double x, double h,
                double xini, const double* __restrict__ H0,
                const double* __restrict__ b1_proj,
                const GPUDensityProfileDevice& profile,
                double HI_constants, bool is_antinu,
                double* __restrict__ y_out)
{
  double HI[9], k[9], acc[9], tmp[9];

  // k1 = f(x, y) = iCommutator(y, HI(x))
  computeHI_SU3(x, xini, H0, b1_proj, profile, HI_constants, is_antinu, HI);
  iCommutatorSU3(y, HI, k);
  #pragma unroll
  for (int c = 0; c < 9; c++) {
    acc[c] = k[c] / 6.0;
    tmp[c] = y[c] + 0.5 * h * k[c];
  }

  // k2 = f(x + h/2, y + h/2*k1) = iCommutator(tmp, HI(x+h/2))
  computeHI_SU3(x + 0.5*h, xini, H0, b1_proj, profile, HI_constants, is_antinu, HI);
  iCommutatorSU3(tmp, HI, k);
  #pragma unroll
  for (int c = 0; c < 9; c++) {
    acc[c] += k[c] / 3.0;
    tmp[c] = y[c] + 0.5 * h * k[c];
  }

  // k3 = f(x + h/2, y + h/2*k2) — reuse HI from k2 (same position)
  iCommutatorSU3(tmp, HI, k);
  #pragma unroll
  for (int c = 0; c < 9; c++) {
    acc[c] += k[c] / 3.0;
    tmp[c] = y[c] + h * k[c];
  }

  // k4 = f(x + h, y + h*k3) = iCommutator(tmp, HI(x+h))
  computeHI_SU3(x + h, xini, H0, b1_proj, profile, HI_constants, is_antinu, HI);
  iCommutatorSU3(tmp, HI, k);
  #pragma unroll
  for (int c = 0; c < 9; c++) {
    acc[c] += k[c] / 6.0;
    y_out[c] = y[c] + h * acc[c];
  }
}

// ============================================================
// Per-substage profile cache (Perf #2).
//
// Layout in shared memory (PROFILE_CACHE_DOUBLES contiguous doubles):
//   [0 .. MAX_TARGETS-1]            : target_ndens[t]   (eV^3, natural units)
//   [MAX_TARGETS .. 2*MAX_TARGETS-1]: target_frac[t]    (= ndens[t] / sum_t ndens)
//   [2*MAX_TARGETS]                 : ye                (electron fraction)
//   [2*MAX_TARGETS + 1]             : num_e             (electron number density, eV^3)
//
// The per-substage time `x_eval` is shared across all threads in the block
// (every thread enters the same RK substage at the same x_sub). Therefore
// each Akima spline only needs to be evaluated *once per substage*, not once
// per thread per call site. Previously target_ndens[16] was recomputed per
// thread inside computeNCCascadeSU3 / computeTauRegenSU3 / computeGlashowCascadeSU3
// / computeDerivativeSU3 — at 128 threads, that's a 128x amplification of the
// same spline work. With substage-level cascade refresh (commit 73eb097)
// this happens 7+ times per macro step, which dominated kernel time.
//
// This cache is populated cooperatively (thread 0 evaluates the splines)
// at the start of each refresh / derivative call site, with a single
// __syncthreads() to publish to the rest of the block.
// ============================================================

static constexpr int PROFILE_CACHE_DOUBLES = 2 * MAX_TARGETS + 2;
static constexpr int CACHE_NDENS_OFFSET = 0;
static constexpr int CACHE_FRAC_OFFSET  = MAX_TARGETS;
static constexpr int CACHE_YE_OFFSET    = 2 * MAX_TARGETS;
static constexpr int CACHE_NUME_OFFSET  = 2 * MAX_TARGETS + 1;

// Populate the per-substage profile cache cooperatively. Caller must follow
// with __syncthreads() before any thread reads from the cache.
// All threads must enter this function (same x_eval, same idata, same profile).
__device__ __forceinline__
void populateProfileCacheSU3(const GPUDensityProfileDevice& profile,
                             double x_eval,
                             const InteractionDataGPU& idata,
                             double* __restrict__ s_cache)
{
  // Thread 0 fills the entire cache. The work is O(n_targets+2) Akima
  // evaluations, dwarfed by the cooperative cascade kernels that follow,
  // so spreading the work across threads adds branching for negligible gain.
  if (threadIdx.x == 0) {
    int n_t = idata.n_targets;
    if (n_t > MAX_TARGETS) n_t = MAX_TARGETS;
    double total = 0.0;
    for (int t = 0; t < n_t; t++) {
      double v = evaluateTargetFraction(profile, t, x_eval);
      s_cache[CACHE_NDENS_OFFSET + t] = v;
      total += v;
    }
    // Zero-fill unused target slots so callers iterating up to n_targets
    // (with n_targets > some prior call's n_targets) can't read stale data.
    for (int t = n_t; t < MAX_TARGETS; t++) {
      s_cache[CACHE_NDENS_OFFSET + t] = 0.0;
      s_cache[CACHE_FRAC_OFFSET + t]  = 0.0;
    }
    double inv = (total > 0.0) ? (1.0 / total) : 0.0;
    for (int t = 0; t < n_t; t++)
      s_cache[CACHE_FRAC_OFFSET + t] = s_cache[CACHE_NDENS_OFFSET + t] * inv;

    double ye = evaluateYe(profile, x_eval);
    s_cache[CACHE_YE_OFFSET] = ye;

    // num_e: matches electronNumberDensitySU3() logic but evaluates once.
    double num_e = 0.0;
    if (profile.has_num_e) {
      num_e = evaluateNumE(profile, x_eval);
    } else if (n_t == 1) {
      num_e = s_cache[CACHE_NDENS_OFFSET + 0] * ye;
    } else if (n_t > 1) {
      num_e = s_cache[CACHE_NDENS_OFFSET + 0]; // proton density (n_e = n_p)
    }
    s_cache[CACHE_NUME_OFFSET] = num_e;
  }
}

// ============================================================
// Cooperative NC cascade computation for SU(3)
//
// Computes nc_factors[rho][alpha][e1] using the state in shared memory.
// Called cooperatively by all threads in a block at the start of each step.
//
// nc_factors[rho][alpha][e1] = sum_{e2>e1} sum_target:
//     targetFrac[t] * Tr(evol_proj[rho][alpha][e2] . state[rho][e2])
//     * invlen_NC[rho][alpha][e2] * delE[e2-1] * dNdE_NC[t][rho][alpha][e2][e1]
//
// s_state: [nrhos][ne][SU] — current state in shared memory
// s_nc_factors: [nrhos][3][ne] — output nc_factors in shared memory
// s_profile_cache: per-substage profile cache (see populateProfileCacheSU3)
// ============================================================

__device__
void computeNCCascadeSU3(int ne, int nrhos, int numneu,
                         double x_eval, double xini,
                         const double* __restrict__ H0_array,
                         const double* __restrict__ b1_proj,
                         const InteractionDataGPU& idata,
                         const GPUDensityProfileDevice& profile,
                         const double* __restrict__ s_state,  // shared: [nrhos][ne][9]
                         double* __restrict__ s_nc_factors,    // shared: [nrhos][3][ne]
                         const double* __restrict__ s_profile_cache)
{
  constexpr int SU = 9;

  // Read precomputed target densities from per-substage profile cache.
  // These were evaluated once at x_eval by populateProfileCacheSU3 and
  // published via __syncthreads() before this function runs.
  const double* target_ndens = s_profile_cache + CACHE_NDENS_OFFSET;
  const double* target_frac  = s_profile_cache + CACHE_FRAC_OFFSET;

  // Each thread handles a subset of output energies e1
  for (int e1 = threadIdx.x; e1 < ne; e1 += blockDim.x) {
    for (int rho = 0; rho < nrhos; rho++) {
      // Evolved projectors at e1's H0 (for flavor extraction at other energies)
      // Actually, for the cascade we need projectors at each e2, not e1.
      // But the projectors in the interaction picture only depend on H0(e) and time.
      // We compute Tr(evol_proj[rho][alpha][e2] . state[rho][e2]) for each e2.

      for (int alpha = 0; alpha < 3; alpha++) {
        double nc_factor = 0.0;

        // Accumulate contributions from higher energies e2 > e1
        for (int e2 = e1 + 1; e2 < ne; e2++) {
          // Get H0 and projectors at energy e2
          const double* H0_e2 = H0_array + (rho * ne + e2) * SU;
          const double* proj_e2 = b1_proj + rho * numneu * SU;

          // Evolve projector for flavor alpha at time x_eval - xini
          double evol_proj_alpha[9];
          {
            const double s3 = 1.7320508075688772;
            double dt = x_eval - xini;
            double w12 = 2.0 * H0_e2[4];
            double w13 = H0_e2[4] + s3 * H0_e2[8];
            double w23 = H0_e2[4] - s3 * H0_e2[8];
            double CX0, SX0, CX1, SX1, CX2, SX2;
            sincos(w12 * dt, &SX0, &CX0);
            sincos(w13 * dt, &SX1, &CX1);
            sincos(w23 * dt, &SX2, &CX2);
            const double* p = proj_e2 + alpha * SU;
            evol_proj_alpha[0] = p[0];
            evol_proj_alpha[1] = CX0*p[1] + SX0*p[3];
            evol_proj_alpha[2] = CX1*p[2] + SX1*p[6];
            evol_proj_alpha[3] = CX0*p[3] - SX0*p[1];
            evol_proj_alpha[4] = p[4];
            evol_proj_alpha[5] = CX2*p[5] - SX2*p[7];
            evol_proj_alpha[6] = CX1*p[6] - SX1*p[2];
            evol_proj_alpha[7] = CX2*p[7] + SX2*p[5];
            evol_proj_alpha[8] = p[8];
          }

          // Flavor flux at e2: Tr(evol_proj[alpha] . state[rho][e2])
          const double* state_e2 = s_state + (rho * ne + e2) * SU;
          double flux = suTrace3(evol_proj_alpha, state_e2);

          // invlen_NC at e2 using proper number densities
          double invlen_NC_e2 = 0.0;
          for (int t = 0; t < idata.n_targets; t++) {
            size_t idx = sigma_index(t, rho, alpha, e2, idata.nrhos, 3, idata.rounded_ne);
            invlen_NC_e2 += target_ndens[t] * idata.d_sigma_NC[idx];
          }

          // Energy bin width
          double dE = idata.d_delE[e2 - 1];

          // Accumulate cascade contribution
          double flux_weighted = flux * invlen_NC_e2 * dE;
          for (int t = 0; t < idata.n_targets; t++) {
            size_t dNdE_idx = dNdE_index(t, rho, alpha, e2, e1,
                                          idata.nrhos, 3, ne, idata.rounded_ne);
            nc_factor += target_frac[t] * flux_weighted * idata.d_dNdE_NC[dNdE_idx];
          }
        }

        s_nc_factors[(rho * 3 + alpha) * ne + e1] = nc_factor;
      }
    }
  }
}

// ============================================================
// Cooperative tau regeneration computation for SU(3)
// Two-pass O(ne²): production of taus from CC, then decay back to neutrinos.
// Adds tau_lep and tau_hadlep contributions to nc_factors in-place.
// ============================================================

__device__
void computeTauRegenSU3(int ne, int nrhos, int numneu,
                        double x_eval, double xini,
                        const double* __restrict__ H0_array,
                        const double* __restrict__ b1_proj,
                        const InteractionDataGPU& idata,
                        const GPUDensityProfileDevice& profile,
                        const double* __restrict__ s_state,
                        double* __restrict__ s_nc_factors,
                        const double* __restrict__ s_profile_cache)
{
  constexpr int SU = 9;
  constexpr int tau_flavor = 2;

  // Read precomputed target densities from per-substage profile cache
  // (matches CPU's GetTargetNumberDensities() in src/nuSQuIDS.cpp:1061-1064).
  // invlen_CC_tau = sum_t target_ndens[t] * sigma_CC[t,rho,tau,en] and
  // targetFractions[t] = target_ndens[t] / total_ndens.
  const double* target_ndens = s_profile_cache + CACHE_NDENS_OFFSET;
  const double* tfrac        = s_profile_cache + CACHE_FRAC_OFFSET;

  // Pass 1: Compute tau decay fluxes — each thread handles a subset of et
  // Use thread-local accumulation, then add to nc_factors atomically
  // tau_decay_fluxes[et] = sum_{en>et} nu_tau_flux(en) * invlen_CC_tau(en) * dE * dNdE_CC(en->et) * dEt
  // For simplicity, we compute this as a single-pass per thread over the (en, e1) space.

  // Each thread handles a subset of final neutrino energies e1
  for (int e1 = threadIdx.x; e1 < ne; e1 += blockDim.x) {
    double tau_lep_nu = 0.0, tau_lep_nubar = 0.0;
    double tau_hadlep_nu = 0.0, tau_hadlep_nubar = 0.0;

    // For each intermediate tau energy et > e1
    for (int et = e1 + 1; et < ne; et++) {
      double tau_flux = 0.0, taubar_flux = 0.0;

      // Accumulate tau production from all en > et
      for (int en = et + 1; en < ne; en++) {
        for (int rho_src = 0; rho_src < nrhos; rho_src++) {
          // Extract tau neutrino flux at en
          const double* H0_en = H0_array + (rho_src * ne + en) * SU;
          const double* proj_en = b1_proj + rho_src * numneu * SU;
          const double* state_en = s_state + (rho_src * ne + en) * SU;

          // Evolve tau projector
          double evol_tau[9];
          {
            const double s3 = 1.7320508075688772;
            double dt = x_eval - xini;
            double w12 = 2.0 * H0_en[4];
            double w13 = H0_en[4] + s3 * H0_en[8];
            double w23 = H0_en[4] - s3 * H0_en[8];
            double CX0, SX0, CX1, SX1, CX2, SX2;
            sincos(w12 * dt, &SX0, &CX0);
            sincos(w13 * dt, &SX1, &CX1);
            sincos(w23 * dt, &SX2, &CX2);
            const double* p = proj_en + tau_flavor * SU;
            evol_tau[0] = p[0]; evol_tau[1] = CX0*p[1] + SX0*p[3];
            evol_tau[2] = CX1*p[2] + SX1*p[6]; evol_tau[3] = CX0*p[3] - SX0*p[1];
            evol_tau[4] = p[4]; evol_tau[5] = CX2*p[5] - SX2*p[7];
            evol_tau[6] = CX1*p[6] - SX1*p[2]; evol_tau[7] = CX2*p[7] + SX2*p[5];
            evol_tau[8] = p[8];
          }

          double flux = suTrace3(evol_tau, state_en);
          // invlen_CC_tau = sum_t target_ndens[t] * sigma_CC[t] in natural
          // units, matching int_state.invlen_CC on CPU.
          double invlen_CC = 0.0;
          for (int t = 0; t < idata.n_targets; t++) {
            size_t idx = sigma_index(t, rho_src, tau_flavor, en, idata.nrhos, 3, idata.rounded_ne);
            invlen_CC += target_ndens[t] * idata.d_sigma_CC[idx];
          }

          double dEn = idata.d_delE[en - 1];
          double dEt = idata.d_delE[et - 1];
          // dNdE_CC[target][rho][tau][en][et]
          for (int t = 0; t < idata.n_targets; t++) {
            size_t dNdE_idx = dNdE_index(t, rho_src, tau_flavor, en, et,
                                          idata.nrhos, 3, ne, idata.rounded_ne);
            double contrib = flux * invlen_CC * dEn * idata.d_dNdE_CC[dNdE_idx] * dEt;
            if (rho_src == 0) tau_flux += tfrac[t] * contrib;
            else              taubar_flux += tfrac[t] * contrib;
          }
        }
      }

      // Decay contributions: tau → neutrinos at energy e1.
      // dNdE_tau_{all,lep}[rho][et][e1] with rounded_ne stride on the fast
      // axis. The rho index labels the source tau species (0 = nu_tau,
      // 1 = nubar_tau). Matching CPU src/nuSQuIDS.cpp:1099-1102:
      //   nu dest, hadlep = nu_tau flux * dNdE_tau_all[0]
      //   nu dest, lep    = nubar_tau flux * dNdE_tau_lep[1]
      //   nubar dest, hadlep = nubar_tau flux * dNdE_tau_all[1]
      //   nubar dest, lep    = nu_tau flux * dNdE_tau_lep[0]
      if (tau_flux > 0.0 || taubar_flux > 0.0) {
        int rne = idata.rounded_ne;
        size_t tau_idx_nu    = (0 * ne + et) * rne + e1; // rho=0 spectrum
        size_t tau_idx_nubar = (1 * ne + et) * rne + e1; // rho=1 spectrum

        tau_hadlep_nu    += tau_flux    * idata.d_dNdE_tau_all[tau_idx_nu];
        tau_lep_nubar    += tau_flux    * idata.d_dNdE_tau_lep[tau_idx_nu];
        tau_hadlep_nubar += taubar_flux * idata.d_dNdE_tau_all[tau_idx_nubar];
        tau_lep_nu       += taubar_flux * idata.d_dNdE_tau_lep[tau_idx_nubar];
      }
    }

    // Add tau regen contributions to nc_factors
    // tau_lep → e and mu flavors, tau_hadlep → tau flavor
    if (nrhos >= 1) {
      s_nc_factors[(0 * 3 + 0) * ne + e1] += tau_lep_nu;    // electron
      s_nc_factors[(0 * 3 + 1) * ne + e1] += tau_lep_nu;    // muon
      s_nc_factors[(0 * 3 + 2) * ne + e1] += tau_hadlep_nu; // tau
    }
    if (nrhos >= 2) {
      s_nc_factors[(1 * 3 + 0) * ne + e1] += tau_lep_nubar;
      s_nc_factors[(1 * 3 + 1) * ne + e1] += tau_lep_nubar;
      s_nc_factors[(1 * 3 + 2) * ne + e1] += tau_hadlep_nubar;
    }
  }
}

// ============================================================
// Cooperative Glashow resonance computation for SU(3)
// Only for electron antineutrinos.
// Adds gr_factors to all flavors of nc_factors.
// ============================================================

__device__
void computeGlashowCascadeSU3(int ne, int nrhos, int numneu,
                               double x_eval, double xini,
                               const double* __restrict__ H0_array,
                               const double* __restrict__ b1_proj,
                               const InteractionDataGPU& idata,
                               const GPUDensityProfileDevice& profile,
                               int NT_type,
                               const double* __restrict__ s_state,
                               double* __restrict__ s_nc_factors,
                               const double* __restrict__ s_profile_cache)
{
  constexpr int SU = 9;
  // Glashow resonance only affects electron antineutrinos
  int rho = (NT_type == 3) ? 1 : 0; // antineutrino rho index

  // Read precomputed electron number density from per-substage profile cache.
  // populateProfileCacheSU3() mirrors electronNumberDensitySU3() (isoscalar,
  // p/n, body-composition, nuclear-XS branches) so this is a drop-in lookup.
  double num_e = s_profile_cache[CACHE_NUME_OFFSET];

  for (int e1 = threadIdx.x; e1 < ne; e1 += blockDim.x) {
    double gr_factor = 0.0;

    for (int e2 = e1 + 1; e2 < ne; e2++) {
      // Extract electron antineutrino flux at e2
      const double* H0_e2 = H0_array + (rho * ne + e2) * SU;
      const double* proj_e2 = b1_proj + rho * numneu * SU;
      const double* state_e2 = s_state + (rho * ne + e2) * SU;

      // Evolve electron projector (flavor 0)
      double evol_e[9];
      {
        const double s3 = 1.7320508075688772;
        double dt = x_eval - xini;
        double w12 = 2.0 * H0_e2[4];
        double w13 = H0_e2[4] + s3 * H0_e2[8];
        double w23 = H0_e2[4] - s3 * H0_e2[8];
        double CX0, SX0, CX1, SX1, CX2, SX2;
        sincos(w12 * dt, &SX0, &CX0);
        sincos(w13 * dt, &SX1, &CX1);
        sincos(w23 * dt, &SX2, &CX2);
        const double* p = proj_e2; // flavor 0 = electron
        evol_e[0] = p[0]; evol_e[1] = CX0*p[1] + SX0*p[3];
        evol_e[2] = CX1*p[2] + SX1*p[6]; evol_e[3] = CX0*p[3] - SX0*p[1];
        evol_e[4] = p[4]; evol_e[5] = CX2*p[5] - SX2*p[7];
        evol_e[6] = CX1*p[6] - SX1*p[2]; evol_e[7] = CX2*p[7] + SX2*p[5];
        evol_e[8] = p[8];
      }

      double flux = suTrace3(evol_e, state_e2);

      // invlen_GR at e2 = sigma_GR * electron number density (natural units).
      double invlen_GR = idata.d_sigma_GR[e2] * num_e;

      double dE = idata.d_delE[e2 - 1];
      gr_factor += flux * invlen_GR * dE * idata.d_dNdE_GR[e2 * idata.rounded_ne + e1];
    }

    // Glashow contributes equally to all flavors (for the antineutrino rho)
    for (int flv = 0; flv < 3; flv++)
      s_nc_factors[(rho * 3 + flv) * ne + e1] += gr_factor;
  }
}

// ============================================================
// Compute InteractionsRho for SU(3)
// InteractionsRho = sum_flv nc_factors[rho][flv][e1] * evol_proj[flv]
// (Tau regeneration and Glashow will be added in Phase 5)
// ============================================================

__device__ __forceinline__
void computeInteractionsRhoSU3(int ie, int rho, int ne,
                                const double* __restrict__ nc_factors, // [nrhos][3][ne]
                                const double* __restrict__ evol_proj,  // [3][9]
                                double* __restrict__ F)                // [9] output
{
  #pragma unroll
  for (int c = 0; c < 9; c++) F[c] = 0.0;

  for (int flv = 0; flv < 3; flv++) {
    double nc_f = nc_factors[(rho * 3 + flv) * ne + ie];
    const double* ep = evol_proj + flv * 9;
    #pragma unroll
    for (int c = 0; c < 9; c++)
      F[c] += nc_f * ep[c];
  }
}

// ============================================================
// Refresh nc_factors / tau_lep / tau_hadlep / Glashow contributions
// using the current `s_state` as the source state at time `x_eval`.
//
// This is cooperative: every thread in the block participates in the
// computeNCCascadeSU3 / computeTauRegenSU3 / computeGlashowCascadeSU3
// helpers (they each loop e1 = threadIdx.x; e1 < ne; e1 += blockDim.x).
//
// The caller is responsible for the __syncthreads() *before* this call
// to publish s_state to all threads; this routine emits the syncs *after*
// each contributor so subsequent reads of s_nc_factors see the updated
// values. Matches the CPU's per-substage PreDerive() refresh.
// ============================================================

__device__ __forceinline__
void refreshCascadeFactorsSU3(int ne, int nrhos, int numneu,
                              double x_eval, double xini,
                              const double* __restrict__ H0_array,
                              const double* __restrict__ b1_proj,
                              const InteractionDataGPU& idata,
                              const GPUDensityProfileDevice& profile,
                              bool tauregeneration,
                              bool iglashow, int NT_type,
                              const double* __restrict__ s_state,
                              double* __restrict__ s_nc_factors,
                              double* __restrict__ s_profile_cache)
{
  // Refresh the per-substage profile cache for this x_eval. All cascade
  // contributors and any subsequent computeDerivativeSU3 call at the same
  // substage time read from this cache instead of recomputing target
  // densities / ye / num_e per thread.
  populateProfileCacheSU3(profile, x_eval, idata, s_profile_cache);
  __syncthreads();

  computeNCCascadeSU3(ne, nrhos, numneu, x_eval, xini,
                      H0_array, b1_proj, idata, profile,
                      s_state, s_nc_factors, s_profile_cache);
  __syncthreads();

  if (tauregeneration && idata.d_dNdE_tau_all) {
    computeTauRegenSU3(ne, nrhos, numneu, x_eval, xini,
                       H0_array, b1_proj, idata, profile,
                       s_state, s_nc_factors, s_profile_cache);
    __syncthreads();
  }

  if (iglashow && idata.d_sigma_GR) {
    computeGlashowCascadeSU3(ne, nrhos, numneu, x_eval, xini,
                              H0_array, b1_proj, idata, profile,
                              NT_type, s_state, s_nc_factors, s_profile_cache);
    __syncthreads();
  }
}

// ============================================================
// Compute full interacting derivative for SU(3)
// dρ/dt = i[ρ, HI] - {Γ, ρ} + F_interactions
//
// Explicit device function — replaces lambda to avoid CUDA
// device lambda capture issues that can silently zero results.
// ============================================================

__device__ __forceinline__
void computeDerivativeSU3(double x_eval, double xini,
                          const double* __restrict__ state,
                          const double* __restrict__ H0,
                          const double* __restrict__ b1_proj,
                          const GPUDensityProfileDevice& profile,
                          double HI_constants, bool is_antinu,
                          int ie, int rho, int ne, int numneu,
                          const InteractionDataGPU& idata,
                          bool iglashow, int NT_type,
                          const double* __restrict__ nc_factors,
                          const double* __restrict__ s_profile_cache,
                          double* __restrict__ deriv)
{
  // Coherent term: reuse the proven computeHI_SU3 + iCommutator path
  double HI[9];
  computeHI_SU3(x_eval, xini, H0, b1_proj, profile, HI_constants, is_antinu, HI);
  iCommutatorSU3(state, HI, deriv);  // deriv = i[ρ, HI]

  // Absorption: -ACommutator(Gamma, ρ)
  // Read target_ndens, ye, num_e from the per-substage profile cache (Perf #2).
  // Caller guarantees the cache was populated for this x_eval.
  const double* target_ndens = s_profile_cache + CACHE_NDENS_OFFSET;
  double ye_x  = s_profile_cache[CACHE_YE_OFFSET];
  double num_e = s_profile_cache[CACHE_NUME_OFFSET];
  double invlen[3];
  computeInvlenSU3(ie, rho, ne, target_ndens, idata,
                   iglashow, NT_type, ye_x, num_e, invlen);

  double evol_proj[3 * 9];
  evolveProjectorsSU3(x_eval, xini, H0, b1_proj, numneu, evol_proj);

  double Gamma[9];
  computeGammaRhoSU3(invlen, evol_proj, Gamma);

  double acomm[9];
  antiCommutatorSU3(Gamma, state, acomm);

  // Cascade source term
  double F_int[9] = {0,0,0,0,0,0,0,0,0};
  if (nc_factors) {
    computeInteractionsRhoSU3(ie, rho, ne, nc_factors, evol_proj, F_int);
  }

  // deriv = i[ρ, HI] - {Γ, ρ} + F_interactions
  #pragma unroll
  for (int c = 0; c < 9; c++)
    deriv[c] += -acomm[c] + F_int[c];
}

// ============================================================
// Main evolve kernel — SU(3)
//
// One block per path (zenith angle).
// Threads parallelize over energy nodes.
// Adaptive step-doubling RK4 in the interaction picture.
//
// Algorithm per step:
//   1. Full RK4 step of size h → sf
//   2. Two half-RK4 steps of size h/2 → sh
//   3. Richardson extrapolation: y_corrected = sh + (sh - sf)/15
//   4. Error estimate: |sf - sh|/15
//   5. PI step-size controller for next step
// ============================================================

template<int NFLV>
__global__
void __launch_bounds__(128)
evolveKernelImpl(const PhysicsParams params,
                 const PathDeviceData* __restrict__ paths,
                 const double* __restrict__ H0_array,
                 const double* __restrict__ b1_proj,
                 const InteractionDataGPU* __restrict__ interaction_data_ptr,
                 double* __restrict__ d_workspace_corrected,
                 double* __restrict__ d_workspace_sf,
                 double* __restrict__ d_workspace_k1,
                 double* __restrict__ d_workspace_k2,
                 double* __restrict__ d_workspace_k3,
                 double* __restrict__ d_workspace_k4,
                 double* __restrict__ d_workspace_k5,
                 double* __restrict__ d_workspace_k6,
                 double* __restrict__ d_workspace_k7,
                 const SolverConfig solver_config,
                 double* __restrict__ states,
                 int n_paths)
{
  // Only SU(3) has a full implementation; SU(4) placeholder below
  if (NFLV != 3) return;

  constexpr int SU = NFLV * NFLV;

  int path_idx = blockIdx.x;
  if (path_idx >= n_paths) return;

  // Copy interaction data to local (avoids repeated global dereference)
  InteractionDataGPU interaction_data;
  if (interaction_data_ptr) {
    interaction_data = *interaction_data_ptr;
  } else {
    interaction_data.n_targets = 0;
    interaction_data.nrhos = 0;
    interaction_data.numneu = 0;
    interaction_data.ne = 0;
    interaction_data.d_dNdE_CC = nullptr;
    interaction_data.d_dNdE_NC = nullptr;
    interaction_data.d_sigma_CC = nullptr;
    interaction_data.d_sigma_NC = nullptr;
    interaction_data.d_dNdE_tau_all = nullptr;
    interaction_data.d_dNdE_tau_lep = nullptr;
    interaction_data.d_sigma_GR = nullptr;
    interaction_data.d_dNdE_GR = nullptr;
    interaction_data.d_energies = nullptr;
    interaction_data.d_delE = nullptr;
    interaction_data.d_b0_proj = nullptr;
    interaction_data.d_b1_proj = nullptr;
    interaction_data.total_bytes = 0;
  }

  const int ne = params.ne;
  const int nrhos = params.nrhos;
  const int state_size = ne * nrhos * SU;

  const PathDeviceData& path = paths[path_idx];
  double* my_state = states + path_idx * state_size;

  double xini = path.xini;
  double xend = path.xend;
  double total_length = xend - xini;

  // Vacuum: in interaction picture, d/dt rho_tilde = 0 → no change
  if (path.profile.type == ProfileType::VACUUM || !params.ioscillations)
    return;
  if (total_length <= 0.0)
    return;

  // Initialize step size proportional to total path length.
  // In the interaction picture, the evolution is driven by the matter
  // potential which varies on density-change scales (~100 km), much slower
  // than vacuum oscillations. Starting at 1% of total path is reasonable.
  double h = total_length * 0.01;
  h = fmax(h, solver_config.h_min);
  h = fmin(h, solver_config.h_max);

  double x = xini;
  int step_count = 0;

  // Per-path slices of the persistent integrator staging workspaces.
  // Layout: [path_idx * nrhos * ne * SU + (rho * ne + ie) * SU + c].
  // This replaces the former per-thread corrected_buf[MAX_PAIRS*SU] +
  // sf_buf[MAX_PAIRS*SU] locals, which were the dominant source of stack
  // spill and register pressure.
  //
  // Oscillation-only path uses path_corrected and path_sf (RK4 + Richardson).
  // Interaction path additionally uses path_k[0..6] for DOPRI5 stage
  // derivatives and path_sf as scratch for the candidate y_{n+1}.
  const int path_state_size = nrhos * ne * SU;
  double* path_corrected = d_workspace_corrected + (size_t)path_idx * path_state_size;
  double* path_sf        = d_workspace_sf        + (size_t)path_idx * path_state_size;
  double* path_k[7];
  if (d_workspace_k1) {
    path_k[0] = d_workspace_k1 + (size_t)path_idx * path_state_size;
    path_k[1] = d_workspace_k2 + (size_t)path_idx * path_state_size;
    path_k[2] = d_workspace_k3 + (size_t)path_idx * path_state_size;
    path_k[3] = d_workspace_k4 + (size_t)path_idx * path_state_size;
    path_k[4] = d_workspace_k5 + (size_t)path_idx * path_state_size;
    path_k[5] = d_workspace_k6 + (size_t)path_idx * path_state_size;
    path_k[6] = d_workspace_k7 + (size_t)path_idx * path_state_size;
  } else {
    #pragma unroll
    for (int i = 0; i < 7; i++) path_k[i] = nullptr;
  }

  // Shared memory layout for interactions (when enabled):
  // [0 .. nrhos*ne*SU)                              : s_state          (cascade source state)
  // [nrhos*ne*SU      .. + nrhos*3*ne)              : s_nc_factors     (cascade output)
  // [+ nrhos*3*ne     .. + PROFILE_CACHE_DOUBLES)   : s_profile_cache  (per-substage profile cache, Perf #2)
  // Dynamic shared memory is allocated in the kernel launch (see launchEvolve).
  extern __shared__ double smem[];
  const bool do_interactions = params.iinteraction && interaction_data.n_targets > 0;
  double* s_state         = smem;                                                    // [nrhos][ne][SU]
  double* s_nc_factors    = smem + nrhos * ne * SU;                                  // [nrhos][3][ne]
  double* s_profile_cache = smem + nrhos * ne * SU + nrhos * 3 * ne;                 // [PROFILE_CACHE_DOUBLES]

  // FSAL (First-Same-As-Last) flag for the DOPRI5 interaction-path integrator.
  // When valid, path_k[0] (== k1 of the upcoming attempted step) already holds
  // f(t_n, y_n) — either as the previous accepted step's k7 or as k1 of the
  // previous (rejected) attempt at the same (t_n, y_n). Invariant: path_k[0]
  // is consistent with the *current* (t_n, y_n) iff fsal_valid == true.
  __shared__ bool fsal_valid;
  if (threadIdx.x == 0) fsal_valid = false;
  __syncthreads();

  while (x < xend - 1.0e-15 * total_length && step_count < solver_config.max_steps) {
    double h_try = fmin(h, xend - x);
    if (h_try <= 0.0) break;

    // ============================================================
    // Cooperative interaction precomputation (Phase 4)
    // Compute nc_factors using the current accepted state.
    // This uses a "lagged" approximation — nc_factors are computed
    // at the start of the step and held constant during sub-steps.
    // The adaptive step controller ensures this approximation is bounded.
    // ============================================================
    if (do_interactions && !fsal_valid) {
      // Load current state to shared memory and refresh cascade factors at
      // time x. This is the source-state for stage 1 (k1) of the upcoming
      // DOPRI5 step.
      //
      // FSAL fast path: when fsal_valid is true, both s_state and
      // s_nc_factors are already consistent with (x, y_n):
      //  - on accept: stage 7 left s_state = y_{n+1} = new y_n, with
      //               cascade factors refreshed at (x_old + h, y_{n+1}) =
      //               (x_new, y_n_new), and path_k[0] holds k1 = stage 7's k7.
      //  - on reject: my_state never changed; we re-run with smaller h, so
      //               stage 1's k = f(x, y_n) is unchanged and path_k[0]
      //               is still valid. The substage chain that follows
      //               will overwrite s_state and s_nc_factors anyway.
      for (int idx = threadIdx.x; idx < nrhos * ne * SU; idx += blockDim.x)
        s_state[idx] = my_state[idx];
      __syncthreads();

      refreshCascadeFactorsSU3(ne, nrhos, NFLV, x, xini,
                                H0_array, b1_proj, interaction_data,
                                path.profile, params.tauregeneration,
                                params.iglashow, params.NT_type,
                                s_state, s_nc_factors, s_profile_cache);
    }

    double local_max_err = 0.0;

    // For oscillation-only, all computation is per-thread with no sync needed.
    // For interactions, we split into phases to refresh nc_factors between
    // the two half-steps, matching the CPU's per-evaluation cascade update.
    if (!do_interactions) {
      // ---- Oscillation-only path: single pass, no sync ----
      for (int ie = threadIdx.x; ie < ne; ie += blockDim.x) {
        for (int rho = 0; rho < nrhos; rho++) {
          double* state_ptr = my_state + (rho * ne + ie) * SU;
          const double* H0 = H0_array + (rho * ne + ie) * SU;
          const double* proj = b1_proj + rho * NFLV * SU;
          bool is_antinu = ((rho == 1) && params.NT_type == 3)
                         || (params.NT_type == 2);

          double y[SU];
          #pragma unroll
          for (int c = 0; c < SU; c++) y[c] = state_ptr[c];

          double sf[SU];
          rk4StepSU3(y, x, h_try, xini, H0, proj, path.profile,
                    params.HI_constants, is_antinu, sf);

          double st[SU], sh[SU];
          rk4StepSU3(y, x, h_try * 0.5, xini, H0, proj, path.profile,
                    params.HI_constants, is_antinu, st);
          rk4StepSU3(st, x + h_try * 0.5, h_try * 0.5, xini, H0, proj, path.profile,
                    params.HI_constants, is_antinu, sh);

          double* corr = path_corrected + (rho * ne + ie) * SU;
          #pragma unroll
          for (int c = 0; c < SU; c++) {
            corr[c] = sh[c] + (sh[c] - sf[c]) / 15.0;
            double err = fabs(sf[c] - sh[c]) / 15.0;
            double scale = solver_config.abs_error +
                           solver_config.rel_error * fmax(fabs(y[c]), fabs(sh[c]));
            if (scale > 0.0)
              local_max_err = fmax(local_max_err, err / scale);
          }
        }
      }
    } else {
      // ---- Interaction path: Dormand-Prince 5(4) embedded pair (DOPRI5).
      //
      //      Replaces the previous step-doubling RK4 + Richardson extrapolation
      //      scheme. DOPRI5 produces a 5th-order solution and a 4th-order
      //      embedded estimate from the same 7 stage derivatives, so we get
      //      both `y_{n+1}` and the truncation-error estimate from a single
      //      pass — no sf/sh redo. The First-Same-As-Last (FSAL) property
      //      lets us reuse stage 7's derivative as the next accepted step's
      //      stage 1, saving one RHS evaluation per accepted step.
      //
      //      The 7 stage derivatives k1..k7 live in per-path slabs of the
      //      Propagator-owned d_workspace_k1..k7 (sized like path_sf), so
      //      they don't bloat per-thread state. Each thread loads k_j[slot]
      //      from global memory when needed and writes its own k_i[slot]
      //      directly to global memory after the substage compute. This
      //      keeps register pressure close to what we had with two RK4
      //      passes (which already needed sf and sh staging).
      //
      //      Cascade refresh count drops from 12 (Richardson) to 7 per
      //      attempted step (1 start-of-step + 5 substage times for stages
      //      2..6 + 1 at endpoint for FSAL stage 7). The endpoint refresh
      //      uses the 5th-order y_{n+1}, so on accept we can reuse k7 as
      //      the next step's k1 without recomputing cascade factors at t_n.
      //
      //      DOPRI5 Butcher tableau:
      //        c[i] : substage time fractions
      //        a[i] : input weights (k1..k_{i} on stage i+1)
      //        b    : 5th-order weights — applied to k1..k6 for y_{n+1}
      //                (k7 weight is 0; k7 only feeds next step's FSAL)
      //        e    : (b - b_hat) — error weights, applied to all 7 k's

      // Per-thread storage. With EVOLVE_THREADS=128 and nrhos<=2 this covers
      // ne up to 256. For ne=200 that's MAX_INT_PAIRS=4 per thread; for
      // typical interaction tests ne<=40, MAX_INT_PAIRS=1 actually suffices
      // but we keep 4 for safety.
      constexpr int MAX_INT_PAIRS = 4;

      int slot_ie[MAX_INT_PAIRS];
      int slot_rho[MAX_INT_PAIRS];
      int n_slots = 0;
      for (int ie = threadIdx.x; ie < ne; ie += blockDim.x) {
        for (int rho = 0; rho < nrhos; rho++) {
          if (n_slots < MAX_INT_PAIRS) {
            slot_ie[n_slots] = ie;
            slot_rho[n_slots] = rho;
          }
          n_slots++;
        }
      }
      if (n_slots > MAX_INT_PAIRS) {
        if (threadIdx.x == 0 && blockIdx.x == 0 && step_count == 0) {
          printf("ERROR: DOPRI5 kernel needs MAX_INT_PAIRS>=%d "
                 "(ne=%d, nrhos=%d, blockDim=%d). Falling back will give "
                 "wrong results.\n", n_slots, ne, nrhos, blockDim.x);
        }
        n_slots = MAX_INT_PAIRS;
      }

      // Per-thread y_init holds the start-of-step state y_n for each slot.
      // Reading y_n[slot] from global my_state once at the top of the macro
      // step is cheaper than reading it 6 times during stage assembly.
      double y_init[MAX_INT_PAIRS][SU];
      for (int s = 0; s < n_slots; s++) {
        const double* state_ptr = my_state + (slot_rho[s] * ne + slot_ie[s]) * SU;
        #pragma unroll
        for (int c = 0; c < SU; c++) y_init[s][c] = state_ptr[c];
      }

      // Dormand-Prince 5(4) Butcher tableau (standard form).
      // c[stage] = substage time as fraction of h_try.
      const double dp_c1 = 1.0/5.0;
      const double dp_c2 = 3.0/10.0;
      const double dp_c3 = 4.0/5.0;
      const double dp_c4 = 8.0/9.0;
      // c5 = c6 = 1.0

      // a[i][j]: weight of k_{j+1} when assembling stage i+1's input state.
      // Only the actually-used coefficients are listed; zeros are skipped.
      // Stage 2 (computes k2): y + h*(1/5)*k1
      const double dp_a21 = 1.0/5.0;
      // Stage 3 (computes k3): y + h*(3/40, 9/40) * (k1, k2)
      const double dp_a31 = 3.0/40.0;
      const double dp_a32 = 9.0/40.0;
      // Stage 4 (computes k4): y + h*(44/45, -56/15, 32/9) * (k1, k2, k3)
      const double dp_a41 = 44.0/45.0;
      const double dp_a42 = -56.0/15.0;
      const double dp_a43 = 32.0/9.0;
      // Stage 5 (computes k5): y + h*(19372/6561, -25360/2187, 64448/6561, -212/729)
      const double dp_a51 = 19372.0/6561.0;
      const double dp_a52 = -25360.0/2187.0;
      const double dp_a53 = 64448.0/6561.0;
      const double dp_a54 = -212.0/729.0;
      // Stage 6 (computes k6): y + h*(9017/3168, -355/33, 46732/5247, 49/176, -5103/18656)
      const double dp_a61 = 9017.0/3168.0;
      const double dp_a62 = -355.0/33.0;
      const double dp_a63 = 46732.0/5247.0;
      const double dp_a64 = 49.0/176.0;
      const double dp_a65 = -5103.0/18656.0;

      // 5th-order solution weights (also = stage-7 input row, since DOPRI5 is
      // FSAL — stage 7 evaluates f at y_{n+1}).
      // y_{n+1} = y_n + h * (b1*k1 + 0*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
      const double dp_b1 = 35.0/384.0;
      const double dp_b3 = 500.0/1113.0;
      const double dp_b4 = 125.0/192.0;
      const double dp_b5 = -2187.0/6784.0;
      const double dp_b6 = 11.0/84.0;

      // Error weights e[i] = b[i] - b_hat[i]. err = h * sum e_i * k_i.
      const double dp_e1 = 71.0/57600.0;
      // dp_e2 = 0.0
      const double dp_e3 = -71.0/16695.0;
      const double dp_e4 = 71.0/1920.0;
      const double dp_e5 = -17253.0/339200.0;
      const double dp_e6 = 22.0/525.0;
      const double dp_e7 = -1.0/40.0;

      // ============================================================
      // Stage 1: compute k1 = f(x, y_n).
      // FSAL: if path_k[0] already holds f(x, y_n) from the previous step
      // (or rejected attempt at the same (x, y_n)), reuse it.
      // ============================================================
      if (!fsal_valid) {
        // s_nc_factors and s_profile_cache were just refreshed at (x, my_state)
        // by the start-of-step block above. Use them to compute k1.
        for (int s = 0; s < n_slots; s++) {
          int ie = slot_ie[s];
          int rho = slot_rho[s];
          const double* H0 = H0_array + (rho * ne + ie) * SU;
          const double* proj = b1_proj + rho * NFLV * SU;
          bool is_antinu = ((rho == 1) && params.NT_type == 3)
                         || (params.NT_type == 2);
          double k_local[SU];
          computeDerivativeSU3(x, xini, y_init[s], H0, proj, path.profile,
                               params.HI_constants, is_antinu,
                               ie, rho, ne, NFLV,
                               interaction_data, params.iglashow, params.NT_type,
                               s_nc_factors, s_profile_cache, k_local);
          double* k1_ptr = path_k[0] + (rho * ne + ie) * SU;
          #pragma unroll
          for (int c = 0; c < SU; c++) k1_ptr[c] = k_local[c];
        }
      }
      // After stage 1, k1 is in path_k[0]. (No global sync needed: each thread
      // wrote only the slabs for slots it owns, and stage 2 only reads each
      // thread's own slots back.)

      // ============================================================
      // Stages 2..6: compute k_i = f(x + c_i*h, y_n + h * sum_{j<i} a_{ij} * k_j).
      // Each stage:
      //   (a) Each thread builds y_eval for its slots from y_init + h*linear
      //       combination of previous k's (read from path_k[j] in global mem).
      //   (b) Each thread writes y_eval into s_state[slot] so that the cooperative
      //       cascade refresh at x_sub sees the substage state.
      //   (c) __syncthreads() — publish s_state.
      //   (d) refreshCascadeFactorsSU3 at x_sub; emits its own __syncthreads().
      //   (e) Each thread evaluates the derivative for its slots and writes
      //       k_i[slot] to path_k[i-1] global memory.
      // ============================================================

      // Stage assembly is repeated across stages 2..6. To avoid CUDA device
      // lambdas (which have historically produced silently wrong results in
      // this kernel — see Phase 4 commentary), the assembly is unrolled by
      // stage index. Each block:
      //   1. Each thread computes y_eval = y_init + h * sum_{j<i} a_ij * k_j
      //      and writes it to s_state[slot] (its own slots only).
      //   2. __syncthreads() to publish s_state.
      //   3. Cooperative refreshCascadeFactorsSU3 at x_sub.
      //   4. Each thread evaluates the RHS for its slots and writes the
      //      per-slot k_i to path_k[stage-1].
      //
      // Internal __syncthreads() inside refreshCascadeFactorsSU3 satisfies
      // the publish-before-read requirement for s_state and s_nc_factors.

#define DOPRI5_EVAL_DERIV_AND_STORE(STAGE_INDEX_1BASED)                         \
  do {                                                                          \
    for (int s = 0; s < n_slots; s++) {                                         \
      int ie = slot_ie[s];                                                      \
      int rho = slot_rho[s];                                                    \
      const double* H0 = H0_array + (rho * ne + ie) * SU;                       \
      const double* proj = b1_proj + rho * NFLV * SU;                           \
      bool is_antinu = ((rho == 1) && params.NT_type == 3)                      \
                     || (params.NT_type == 2);                                  \
      double y_eval[SU];                                                        \
      const double* s_ptr = s_state + (rho * ne + ie) * SU;                     \
      _Pragma("unroll")                                                         \
      for (int c = 0; c < SU; c++) y_eval[c] = s_ptr[c];                        \
      double k_local[SU];                                                       \
      computeDerivativeSU3(x_sub, xini, y_eval, H0, proj, path.profile,         \
                           params.HI_constants, is_antinu,                      \
                           ie, rho, ne, NFLV,                                   \
                           interaction_data, params.iglashow, params.NT_type,   \
                           s_nc_factors, s_profile_cache, k_local);             \
      double* k_ptr = path_k[(STAGE_INDEX_1BASED) - 1] + (rho * ne + ie) * SU;  \
      _Pragma("unroll")                                                         \
      for (int c = 0; c < SU; c++) k_ptr[c] = k_local[c];                       \
    }                                                                           \
  } while (0)

      // Stage 2: y_eval = y_init + h * a21 * k1
      {
        const double x_sub = x + dp_c1 * h_try;
        for (int s = 0; s < n_slots; s++) {
          int ie = slot_ie[s], rho = slot_rho[s];
          const double* k1_ptr = path_k[0] + (rho * ne + ie) * SU;
          double* s_ptr = s_state + (rho * ne + ie) * SU;
          #pragma unroll
          for (int c = 0; c < SU; c++)
            s_ptr[c] = y_init[s][c] + h_try * (dp_a21 * k1_ptr[c]);
        }
        __syncthreads();
        refreshCascadeFactorsSU3(ne, nrhos, NFLV, x_sub, xini,
                                  H0_array, b1_proj, interaction_data,
                                  path.profile, params.tauregeneration,
                                  params.iglashow, params.NT_type,
                                  s_state, s_nc_factors, s_profile_cache);
        DOPRI5_EVAL_DERIV_AND_STORE(2);
      }

      // Stage 3: y_eval = y_init + h * (a31*k1 + a32*k2)
      {
        const double x_sub = x + dp_c2 * h_try;
        for (int s = 0; s < n_slots; s++) {
          int ie = slot_ie[s], rho = slot_rho[s];
          const double* k1_ptr = path_k[0] + (rho * ne + ie) * SU;
          const double* k2_ptr = path_k[1] + (rho * ne + ie) * SU;
          double* s_ptr = s_state + (rho * ne + ie) * SU;
          #pragma unroll
          for (int c = 0; c < SU; c++)
            s_ptr[c] = y_init[s][c] + h_try * (dp_a31 * k1_ptr[c]
                                              + dp_a32 * k2_ptr[c]);
        }
        __syncthreads();
        refreshCascadeFactorsSU3(ne, nrhos, NFLV, x_sub, xini,
                                  H0_array, b1_proj, interaction_data,
                                  path.profile, params.tauregeneration,
                                  params.iglashow, params.NT_type,
                                  s_state, s_nc_factors, s_profile_cache);
        DOPRI5_EVAL_DERIV_AND_STORE(3);
      }

      // Stage 4: y_eval = y_init + h * (a41*k1 + a42*k2 + a43*k3)
      {
        const double x_sub = x + dp_c3 * h_try;
        for (int s = 0; s < n_slots; s++) {
          int ie = slot_ie[s], rho = slot_rho[s];
          const double* k1_ptr = path_k[0] + (rho * ne + ie) * SU;
          const double* k2_ptr = path_k[1] + (rho * ne + ie) * SU;
          const double* k3_ptr = path_k[2] + (rho * ne + ie) * SU;
          double* s_ptr = s_state + (rho * ne + ie) * SU;
          #pragma unroll
          for (int c = 0; c < SU; c++)
            s_ptr[c] = y_init[s][c] + h_try * (dp_a41 * k1_ptr[c]
                                              + dp_a42 * k2_ptr[c]
                                              + dp_a43 * k3_ptr[c]);
        }
        __syncthreads();
        refreshCascadeFactorsSU3(ne, nrhos, NFLV, x_sub, xini,
                                  H0_array, b1_proj, interaction_data,
                                  path.profile, params.tauregeneration,
                                  params.iglashow, params.NT_type,
                                  s_state, s_nc_factors, s_profile_cache);
        DOPRI5_EVAL_DERIV_AND_STORE(4);
      }

      // Stage 5: y_eval = y_init + h * (a51*k1 + a52*k2 + a53*k3 + a54*k4)
      {
        const double x_sub = x + dp_c4 * h_try;
        for (int s = 0; s < n_slots; s++) {
          int ie = slot_ie[s], rho = slot_rho[s];
          const double* k1_ptr = path_k[0] + (rho * ne + ie) * SU;
          const double* k2_ptr = path_k[1] + (rho * ne + ie) * SU;
          const double* k3_ptr = path_k[2] + (rho * ne + ie) * SU;
          const double* k4_ptr = path_k[3] + (rho * ne + ie) * SU;
          double* s_ptr = s_state + (rho * ne + ie) * SU;
          #pragma unroll
          for (int c = 0; c < SU; c++)
            s_ptr[c] = y_init[s][c] + h_try * (dp_a51 * k1_ptr[c]
                                              + dp_a52 * k2_ptr[c]
                                              + dp_a53 * k3_ptr[c]
                                              + dp_a54 * k4_ptr[c]);
        }
        __syncthreads();
        refreshCascadeFactorsSU3(ne, nrhos, NFLV, x_sub, xini,
                                  H0_array, b1_proj, interaction_data,
                                  path.profile, params.tauregeneration,
                                  params.iglashow, params.NT_type,
                                  s_state, s_nc_factors, s_profile_cache);
        DOPRI5_EVAL_DERIV_AND_STORE(5);
      }

      // Stage 6: y_eval = y_init + h * (a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5)
      {
        const double x_sub = x + h_try;  // c5 = 1
        for (int s = 0; s < n_slots; s++) {
          int ie = slot_ie[s], rho = slot_rho[s];
          const double* k1_ptr = path_k[0] + (rho * ne + ie) * SU;
          const double* k2_ptr = path_k[1] + (rho * ne + ie) * SU;
          const double* k3_ptr = path_k[2] + (rho * ne + ie) * SU;
          const double* k4_ptr = path_k[3] + (rho * ne + ie) * SU;
          const double* k5_ptr = path_k[4] + (rho * ne + ie) * SU;
          double* s_ptr = s_state + (rho * ne + ie) * SU;
          #pragma unroll
          for (int c = 0; c < SU; c++)
            s_ptr[c] = y_init[s][c] + h_try * (dp_a61 * k1_ptr[c]
                                              + dp_a62 * k2_ptr[c]
                                              + dp_a63 * k3_ptr[c]
                                              + dp_a64 * k4_ptr[c]
                                              + dp_a65 * k5_ptr[c]);
        }
        __syncthreads();
        refreshCascadeFactorsSU3(ne, nrhos, NFLV, x_sub, xini,
                                  H0_array, b1_proj, interaction_data,
                                  path.profile, params.tauregeneration,
                                  params.iglashow, params.NT_type,
                                  s_state, s_nc_factors, s_profile_cache);
        DOPRI5_EVAL_DERIV_AND_STORE(6);
      }

#undef DOPRI5_EVAL_DERIV_AND_STORE

      // ============================================================
      // Stage 7 (FSAL): k7 = f(x + h, y_{n+1}).
      // Combines acceptance candidate y_{n+1} (5th-order) with the same
      // cooperative-refresh chain as stages 2..6. y_{n+1} differs from
      // stage-6 input (which used a[6] = b coincidentally is *not* true
      // in general), so we must redo the cascade refresh at x+h with the
      // accepted y_{n+1} as the source state. We also write y_{n+1} to
      // path_corrected here so the accept branch can simply copy it back.
      // ============================================================
      for (int s = 0; s < n_slots; s++) {
        int ie = slot_ie[s];
        int rho = slot_rho[s];
        const double* k1_ptr = path_k[0] + (rho * ne + ie) * SU;
        const double* k3_ptr = path_k[2] + (rho * ne + ie) * SU;
        const double* k4_ptr = path_k[3] + (rho * ne + ie) * SU;
        const double* k5_ptr = path_k[4] + (rho * ne + ie) * SU;
        const double* k6_ptr = path_k[5] + (rho * ne + ie) * SU;
        double* corr_ptr = path_corrected + (rho * ne + ie) * SU;
        double* s_ptr    = s_state        + (rho * ne + ie) * SU;
        #pragma unroll
        for (int c = 0; c < SU; c++) {
          double y_next = y_init[s][c] + h_try *
            (dp_b1 * k1_ptr[c] + dp_b3 * k3_ptr[c] + dp_b4 * k4_ptr[c]
             + dp_b5 * k5_ptr[c] + dp_b6 * k6_ptr[c]);
          corr_ptr[c] = y_next;
          s_ptr[c]    = y_next;
        }
      }
      __syncthreads();

      refreshCascadeFactorsSU3(ne, nrhos, NFLV, x + h_try, xini,
                                H0_array, b1_proj, interaction_data,
                                path.profile, params.tauregeneration,
                                params.iglashow, params.NT_type,
                                s_state, s_nc_factors, s_profile_cache);

      for (int s = 0; s < n_slots; s++) {
        int ie = slot_ie[s];
        int rho = slot_rho[s];
        const double* H0 = H0_array + (rho * ne + ie) * SU;
        const double* proj = b1_proj + rho * NFLV * SU;
        bool is_antinu = ((rho == 1) && params.NT_type == 3)
                       || (params.NT_type == 2);
        double y_eval[SU];
        const double* s_ptr = s_state + (rho * ne + ie) * SU;
        #pragma unroll
        for (int c = 0; c < SU; c++) y_eval[c] = s_ptr[c];
        double k_local[SU];
        computeDerivativeSU3(x + h_try, xini, y_eval, H0, proj, path.profile,
                             params.HI_constants, is_antinu,
                             ie, rho, ne, NFLV,
                             interaction_data, params.iglashow, params.NT_type,
                             s_nc_factors, s_profile_cache, k_local);
        double* k7_ptr = path_k[6] + (rho * ne + ie) * SU;
        #pragma unroll
        for (int c = 0; c < SU; c++) k7_ptr[c] = k_local[c];
      }

      // ============================================================
      // Embedded error estimate.
      //   err[c] = h * (e1*k1 + e3*k3 + e4*k4 + e5*k5 + e6*k6 + e7*k7)
      //   err_norm = max over (slot, c) of |err[c]| / scale[c]
      //   scale[c] = abs_tol + rel_tol * max(|y_n[c]|, |y_{n+1}[c]|)
      // ============================================================
      for (int s = 0; s < n_slots; s++) {
        int ie = slot_ie[s];
        int rho = slot_rho[s];
        const double* k1_ptr = path_k[0] + (rho * ne + ie) * SU;
        const double* k3_ptr = path_k[2] + (rho * ne + ie) * SU;
        const double* k4_ptr = path_k[3] + (rho * ne + ie) * SU;
        const double* k5_ptr = path_k[4] + (rho * ne + ie) * SU;
        const double* k6_ptr = path_k[5] + (rho * ne + ie) * SU;
        const double* k7_ptr = path_k[6] + (rho * ne + ie) * SU;
        const double* corr_ptr = path_corrected + (rho * ne + ie) * SU;
        #pragma unroll
        for (int c = 0; c < SU; c++) {
          double err = h_try * (dp_e1 * k1_ptr[c] + dp_e3 * k3_ptr[c]
                              + dp_e4 * k4_ptr[c] + dp_e5 * k5_ptr[c]
                              + dp_e6 * k6_ptr[c] + dp_e7 * k7_ptr[c]);
          double y_n   = y_init[s][c];
          double y_np1 = corr_ptr[c];
          double scale = solver_config.abs_error +
                         solver_config.rel_error *
                         fmax(fabs(y_n), fabs(y_np1));
          if (scale > 0.0)
            local_max_err = fmax(local_max_err, fabs(err) / scale);
        }
      }
    }

    // Phase 2: Block-wide error reduction
    double max_err = blockReduceMax(local_max_err);

    // Decide accept/reject via shared memory
    __shared__ bool step_accepted;
    if (threadIdx.x == 0)
      step_accepted = (max_err <= 1.0);
    __syncthreads();

    if (step_accepted) {
      // Phase 3: Write corrected states to global memory.
      for (int ie = threadIdx.x; ie < ne; ie += blockDim.x) {
        for (int rho = 0; rho < nrhos; rho++) {
          double* state_ptr = my_state + (rho * ne + ie) * SU;
          const double* corr = path_corrected + (rho * ne + ie) * SU;
          #pragma unroll
          for (int c = 0; c < SU; c++)
            state_ptr[c] = corr[c];
        }
      }
      x += h_try;
      step_count++;

      // FSAL bookkeeping (interaction path only): the accepted step's k7
      // is f(x_new, y_new), which is exactly what the next step's k1
      // requires. Copy k7 -> k1 cooperatively so subsequent stage 2..6
      // computes can read k1 from path_k[0].
      if (do_interactions) {
        for (int idx = threadIdx.x; idx < nrhos * ne * SU; idx += blockDim.x) {
          path_k[0][idx] = path_k[6][idx];
        }
        if (threadIdx.x == 0) fsal_valid = true;
      }
    } else if (do_interactions) {
      // Step rejected: y_n is unchanged, so path_k[0] (which equals f(x, y_n)
      // computed during this attempt) remains valid for the next attempt at
      // smaller h. k1 only depends on (x, y_n), not on h.
      if (threadIdx.x == 0) fsal_valid = true;
    }
    __syncthreads();

    // PI step-size controller. DOPRI5 5(4) error estimate is 4th-order; the
    // accepted output is 5th-order, so the standard step-size exponent is
    // -1/p with p=5 (the accepted order). Fac=0.9, facmin=0.2, facmax=5.0
    // per the prompt's prescription.
    double factor;
    if (max_err > 1.0e-10)
      factor = 0.9 * pow(max_err, -0.2);
    else
      factor = 5.0;
    factor = fmax(0.2, fmin(5.0, factor));

    h = h_try * factor;
    h = fmax(h, solver_config.h_min);
    h = fmin(h, solver_config.h_max);
  }

  // Debug: report if max steps reached
  if (path_idx == 0 && threadIdx.x == 0) {
    if (step_count >= solver_config.max_steps)
      printf("WARNING: max_steps=%d reached at x=%e (xend=%e) h=%e do_int=%d\n",
             solver_config.max_steps, x, xend, h, (int)do_interactions);
    else
      printf("OK: completed in %d steps, x=%e xend=%e\n", step_count, x, xend);
  }
}

// ============================================================
// Kernel launch wrappers (dispatching by numneu)
// ============================================================

void launchInitH0(const double* d_energies, const double* d_dm2,
                  int ne, int nrhos, int numneu,
                  double* d_H0_array, cudaStream_t stream) {
  // H0 is pre-computed on CPU and uploaded; this kernel is unused
  // but kept for potential future use
  NUSQUIDS_CUDA_CHECK(cudaGetLastError());
}

void launchSetInitialStates(const double* d_flux_data,
                            const double* d_b0_proj,
                            int ne, int nrhos, int numneu,
                            double* d_states,
                            int n_paths, cudaStream_t stream) {
  int su_size = numneu * numneu;
  size_t bytes = (size_t)n_paths * ne * nrhos * su_size * sizeof(double);
  NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(d_states, d_flux_data, bytes,
                                       cudaMemcpyHostToDevice, stream));
}

void launchEvolve(const PhysicsParams& params,
                  const PathDeviceData* d_paths,
                  const double* d_H0_array,
                  const double* d_b1_proj,
                  const InteractionDataGPU* d_interaction_data,
                  double* d_workspace_corrected,
                  double* d_workspace_sf,
                  double* const* d_workspace_k,
                  const SolverConfig& solver_config,
                  double* d_states,
                  int n_paths, int numneu,
                  cudaStream_t stream) {
  if (n_paths == 0) return;

  constexpr int threads = EVOLVE_THREADS;

  // Compute shared memory for interaction cascade computation
  // Layout: s_state[nrhos*ne*SU] + s_nc_factors[nrhos*3*ne] + s_profile_cache[PROFILE_CACHE_DOUBLES]
  // The s_profile_cache (Perf #2) holds target_ndens, target_frac, ye, num_e
  // evaluated once per substage time and broadcast to all threads.
  const bool do_interactions =
      params.iinteraction && d_interaction_data && d_interaction_data->n_targets > 0;
  size_t shared_bytes = 0;
  if (do_interactions) {
    int su_size = numneu * numneu;
    shared_bytes = (size_t)(params.nrhos * params.ne * su_size   // state
                          + params.nrhos * 3 * params.ne         // nc_factors
                          + PROFILE_CACHE_DOUBLES)               // profile cache
                   * sizeof(double);
  }

  // Opt in to the larger per-block shared-memory carveout when the dynamic
  // shared bytes exceed the default 48 KB. On A100 the maximum is 164 KB.
  // This is a no-op on devices that already map the static reservation; the
  // request only matters when shared_bytes > 48*1024.
  static constexpr int A100_MAX_SHARED_BYTES = 165536; // 164 KB - 32 B for system reservation
  if (shared_bytes > 48 * 1024) {
    cudaFuncSetAttribute((const void*)&evolveKernelImpl<3>,
                         cudaFuncAttributeMaxDynamicSharedMemorySize,
                         A100_MAX_SHARED_BYTES);
    cudaFuncSetAttribute((const void*)&evolveKernelImpl<4>,
                         cudaFuncAttributeMaxDynamicSharedMemorySize,
                         A100_MAX_SHARED_BYTES);
    // Clear any error if the device doesn't support this attribute (e.g. MIG).
    cudaGetLastError();
  }

  // Upload InteractionDataGPU struct to device if needed
  // (the struct itself contains device pointers, but it needs to be device-resident
  //  for the kernel to read it via a device pointer)
  InteractionDataGPU* d_idata_on_device = nullptr;
  if (do_interactions) {
    NUSQUIDS_CUDA_CHECK(cudaMalloc(&d_idata_on_device, sizeof(InteractionDataGPU)));
    NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(d_idata_on_device, d_interaction_data,
                                         sizeof(InteractionDataGPU),
                                         cudaMemcpyHostToDevice, stream));
  }

  // DOPRI5 stage-derivative slabs. When interactions are disabled, the kernel
  // never accesses these (do_interactions==false branch is RK4+Richardson),
  // so passing nullptr is safe.
  double* k1 = do_interactions ? d_workspace_k[0] : nullptr;
  double* k2 = do_interactions ? d_workspace_k[1] : nullptr;
  double* k3 = do_interactions ? d_workspace_k[2] : nullptr;
  double* k4 = do_interactions ? d_workspace_k[3] : nullptr;
  double* k5 = do_interactions ? d_workspace_k[4] : nullptr;
  double* k6 = do_interactions ? d_workspace_k[5] : nullptr;
  double* k7 = do_interactions ? d_workspace_k[6] : nullptr;

  // Clear any stale CUDA errors (e.g., from cudaDeviceSetLimit on MIG)
  cudaGetLastError();

  switch (numneu) {
    case 3:
      evolveKernelImpl<3><<<n_paths, threads, shared_bytes, stream>>>(
        params, d_paths, d_H0_array, d_b1_proj, d_idata_on_device,
        d_workspace_corrected, d_workspace_sf,
        k1, k2, k3, k4, k5, k6, k7,
        solver_config, d_states, n_paths);
      break;
    case 4:
      evolveKernelImpl<4><<<n_paths, threads, shared_bytes, stream>>>(
        params, d_paths, d_H0_array, d_b1_proj, d_idata_on_device,
        d_workspace_corrected, d_workspace_sf,
        k1, k2, k3, k4, k5, k6, k7,
        solver_config, d_states, n_paths);
      break;
  }
  NUSQUIDS_CUDA_CHECK(cudaGetLastError());

  // Free the temporary device struct
  if (d_idata_on_device) {
    cudaStreamSynchronize(stream);
    cudaFree(d_idata_on_device);
  }
}

void launchEvalFlavors(const double* d_states,
                       const double* d_evol_proj,
                       int ne, int nrhos, int numneu,
                       double* d_flavors,
                       int n_paths, cudaStream_t stream) {
  // Flavor extraction is done on CPU after state download
}

}} // namespace nusquids::cuda
