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
void __launch_bounds__(128, 2)
evolveKernelImpl(const PhysicsParams params,
                 const PathDeviceData* __restrict__ paths,
                 const double* __restrict__ H0_array,
                 const double* __restrict__ b1_proj,
                 const double* __restrict__ interaction_data,
                 const SolverConfig solver_config,
                 double* __restrict__ states,
                 int n_paths)
{
  // Only SU(3) has a full implementation; SU(4) placeholder below
  if (NFLV != 3) return;

  constexpr int SU = NFLV * NFLV;

  int path_idx = blockIdx.x;
  if (path_idx >= n_paths) return;

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

  // Buffer for corrected states: each thread handles at most
  // ceil(ne/blockDim) * nrhos pairs, each with SU components.
  // Max buffer: 32 pairs * 9 components = 288 doubles.
  constexpr int MAX_PAIRS = 32;
  double corrected_buf[MAX_PAIRS * SU];

  while (x < xend - 1.0e-15 * total_length && step_count < solver_config.max_steps) {
    double h_try = fmin(h, xend - x);
    if (h_try <= 0.0) break;

    double local_max_err = 0.0;
    int pair_idx = 0;

    // Phase 1: Compute corrected states and error for all (ie, rho) pairs.
    // Store corrected states in local buffer — do NOT write to global memory yet.
    for (int ie = threadIdx.x; ie < ne; ie += blockDim.x) {
      for (int rho = 0; rho < nrhos; rho++) {
        double* state_ptr = my_state + (rho * ne + ie) * SU;
        const double* H0 = H0_array + (rho * ne + ie) * SU;
        const double* proj = b1_proj + rho * NFLV * SU;

        bool is_antinu = ((rho == 1) && params.NT_type == 3)
                       || (params.NT_type == 2);

        // Load current state from global memory
        double y[SU];
        #pragma unroll
        for (int c = 0; c < SU; c++) y[c] = state_ptr[c];

        // Full step of size h_try
        double sf[SU];
        rk4StepSU3(y, x, h_try, xini, H0, proj, path.profile,
                    params.HI_constants, is_antinu, sf);

        // Two half-steps of size h_try/2
        double st[SU], sh[SU];
        rk4StepSU3(y, x, h_try * 0.5, xini, H0, proj, path.profile,
                    params.HI_constants, is_antinu, st);
        rk4StepSU3(st, x + h_try * 0.5, h_try * 0.5, xini, H0, proj, path.profile,
                    params.HI_constants, is_antinu, sh);

        // Richardson extrapolation and error estimate
        double* corr = corrected_buf + pair_idx * SU;
        #pragma unroll
        for (int c = 0; c < SU; c++) {
          corr[c] = sh[c] + (sh[c] - sf[c]) / 15.0;
          double err = fabs(sf[c] - sh[c]) / 15.0;
          double scale = solver_config.abs_error +
                         solver_config.rel_error * fmax(fabs(y[c]), fabs(sh[c]));
          if (scale > 0.0)
            local_max_err = fmax(local_max_err, err / scale);
        }
        pair_idx++;
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
      // Phase 3: Write corrected states to global memory
      pair_idx = 0;
      for (int ie = threadIdx.x; ie < ne; ie += blockDim.x) {
        for (int rho = 0; rho < nrhos; rho++) {
          double* state_ptr = my_state + (rho * ne + ie) * SU;
          const double* corr = corrected_buf + pair_idx * SU;
          #pragma unroll
          for (int c = 0; c < SU; c++)
            state_ptr[c] = corr[c];
          pair_idx++;
        }
      }
      x += h_try;
      step_count++;
    }

    // PI step-size controller
    double factor;
    if (max_err > 1.0e-10)
      factor = 0.9 * pow(max_err, -0.2);
    else
      factor = 5.0;
    factor = fmax(0.1, fmin(5.0, factor));

    h = h_try * factor;
    h = fmax(h, solver_config.h_min);
    h = fmin(h, solver_config.h_max);
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
                  const double* d_interaction_data,
                  const SolverConfig& solver_config,
                  double* d_states,
                  int n_paths, int numneu,
                  cudaStream_t stream) {
  if (n_paths == 0) return;

  int threads = 128;

  switch (numneu) {
    case 3:
      evolveKernelImpl<3><<<n_paths, threads, 0, stream>>>(
        params, d_paths, d_H0_array, d_b1_proj,
        d_interaction_data, solver_config, d_states, n_paths);
      break;
    case 4:
      evolveKernelImpl<4><<<n_paths, threads, 0, stream>>>(
        params, d_paths, d_H0_array, d_b1_proj,
        d_interaction_data, solver_config, d_states, n_paths);
      break;
  }
  NUSQUIDS_CUDA_CHECK(cudaGetLastError());
}

void launchEvalFlavors(const double* d_states,
                       const double* d_evol_proj,
                       int ne, int nrhos, int numneu,
                       double* d_flavors,
                       int n_paths, cudaStream_t stream) {
  // Flavor extraction is done on CPU after state download
}

}} // namespace nusquids::cuda
