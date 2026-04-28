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
                      const GPUDensityProfileDevice& profile,
                      double x_eval,
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
      double num_e = electronNumberDensitySU3(profile, x_eval,
                                              target_ndens, idata.n_targets, ye);
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
// ============================================================

__device__
void computeNCCascadeSU3(int ne, int nrhos, int numneu,
                         double x_eval, double xini,
                         const double* __restrict__ H0_array,
                         const double* __restrict__ b1_proj,
                         const InteractionDataGPU& idata,
                         const GPUDensityProfileDevice& profile,
                         const double* __restrict__ s_state,  // shared: [nrhos][ne][9]
                         double* __restrict__ s_nc_factors)    // shared: [nrhos][3][ne]
{
  constexpr int SU = 9;

  // Evaluate precomputed target number densities from profile splines
  // These were computed on the CPU using GetTargetNumberDensities() with squids::Const
  double target_ndens[MAX_TARGETS];
  double target_frac[MAX_TARGETS];
  double total_ndens = 0.0;
  for (int t = 0; t < idata.n_targets && t < MAX_TARGETS; t++) {
    target_ndens[t] = evaluateTargetFraction(profile, t, x_eval);
    total_ndens += target_ndens[t];
  }
  for (int t = 0; t < idata.n_targets && t < MAX_TARGETS; t++)
    target_frac[t] = (total_ndens > 0) ? target_ndens[t] / total_ndens : 0.0;

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
                        double* __restrict__ s_nc_factors)
{
  constexpr int SU = 9;
  constexpr int tau_flavor = 2;

  // Evaluate precomputed target number densities from profile splines
  // (same as computeNCCascadeSU3, matching CPU's GetTargetNumberDensities()).
  // invlen_CC_tau = sum_t target_ndens[t] * sigma_CC[t,rho,tau,en] and
  // targetFractions[t] = target_ndens[t] / total_ndens, matching the CPU
  // code in src/nuSQuIDS.cpp:1061-1064.
  double target_ndens[MAX_TARGETS];
  double tfrac[MAX_TARGETS];
  double total_ndens = 0.0;
  for (int t = 0; t < idata.n_targets && t < MAX_TARGETS; t++) {
    target_ndens[t] = evaluateTargetFraction(profile, t, x_eval);
    total_ndens += target_ndens[t];
  }
  for (int t = 0; t < idata.n_targets && t < MAX_TARGETS; t++)
    tfrac[t] = (total_ndens > 0.0) ? target_ndens[t] / total_ndens : 0.0;

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
                               double* __restrict__ s_nc_factors)
{
  constexpr int SU = 9;
  // Glashow resonance only affects electron antineutrinos
  int rho = (NT_type == 3) ? 1 : 0; // antineutrino rho index

  // Use precomputed electron number density (CPU side mirrors UpdateInteractions
  // and handles isoscalar, p/n, body-composition, and nuclear-XS branches in
  // squids::Const natural units). Falls back to the legacy isoscalar / p-n
  // derivation when num_e is not provided.
  double target_ndens[MAX_TARGETS];
  for (int t = 0; t < idata.n_targets && t < MAX_TARGETS; t++)
    target_ndens[t] = evaluateTargetFraction(profile, t, x_eval);
  double ye = evaluateYe(profile, x_eval);
  double num_e = electronNumberDensitySU3(profile, x_eval,
                                          target_ndens, idata.n_targets, ye);

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
                              double* __restrict__ s_nc_factors)
{
  computeNCCascadeSU3(ne, nrhos, numneu, x_eval, xini,
                      H0_array, b1_proj, idata, profile,
                      s_state, s_nc_factors);
  __syncthreads();

  if (tauregeneration && idata.d_dNdE_tau_all) {
    computeTauRegenSU3(ne, nrhos, numneu, x_eval, xini,
                       H0_array, b1_proj, idata, profile,
                       s_state, s_nc_factors);
    __syncthreads();
  }

  if (iglashow && idata.d_sigma_GR) {
    computeGlashowCascadeSU3(ne, nrhos, numneu, x_eval, xini,
                              H0_array, b1_proj, idata, profile,
                              NT_type, s_state, s_nc_factors);
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
                          double* __restrict__ deriv)
{
  // Coherent term: reuse the proven computeHI_SU3 + iCommutator path
  double HI[9];
  computeHI_SU3(x_eval, xini, H0, b1_proj, profile, HI_constants, is_antinu, HI);
  iCommutatorSU3(state, HI, deriv);  // deriv = i[ρ, HI]

  // Absorption: -ACommutator(Gamma, ρ)
  // Evaluate precomputed target number densities from profile splines
  double target_ndens[MAX_TARGETS];
  for (int t = 0; t < idata.n_targets && t < MAX_TARGETS; t++)
    target_ndens[t] = evaluateTargetFraction(profile, t, x_eval);

  double ye_x = evaluateYe(profile, x_eval);
  double invlen[3];
  computeInvlenSU3(ie, rho, ne, target_ndens, idata,
                   iglashow, NT_type, ye_x, profile, x_eval, invlen);

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
// RK4 step for SU(3) with full interactions
//
// dρ/dt = i[ρ, HI] - {Γ, ρ} + F_interactions
//
// Uses precomputed nc_factors (lagged from start of step).
// ============================================================

__device__
void rk4StepSU3_interacting(const double* __restrict__ y, double x, double h,
                             double xini,
                             const double* __restrict__ H0,
                             const double* __restrict__ b1_proj,
                             const GPUDensityProfileDevice& profile,
                             double HI_constants, bool is_antinu,
                             int ie, int rho, int ne, int numneu,
                             const InteractionDataGPU& idata,
                             bool iglashow, int NT_type,
                             const double* __restrict__ nc_factors,
                             double* __restrict__ y_out)
{
  double k[9], acc[9], tmp[9];

  // k1 = f(x, y)
  computeDerivativeSU3(x, xini, y, H0, b1_proj, profile,
                       HI_constants, is_antinu, ie, rho, ne, numneu,
                       idata, iglashow, NT_type, nc_factors, k);
  #pragma unroll
  for (int c = 0; c < 9; c++) {
    acc[c] = k[c] / 6.0;
    tmp[c] = y[c] + 0.5 * h * k[c];
  }

  // k2 = f(x + h/2, y + h/2*k1)
  computeDerivativeSU3(x + 0.5*h, xini, tmp, H0, b1_proj, profile,
                       HI_constants, is_antinu, ie, rho, ne, numneu,
                       idata, iglashow, NT_type, nc_factors, k);
  #pragma unroll
  for (int c = 0; c < 9; c++) {
    acc[c] += k[c] / 3.0;
    tmp[c] = y[c] + 0.5 * h * k[c];
  }

  // k3 = f(x + h/2, y + h/2*k2)
  computeDerivativeSU3(x + 0.5*h, xini, tmp, H0, b1_proj, profile,
                       HI_constants, is_antinu, ie, rho, ne, numneu,
                       idata, iglashow, NT_type, nc_factors, k);
  #pragma unroll
  for (int c = 0; c < 9; c++) {
    acc[c] += k[c] / 3.0;
    tmp[c] = y[c] + h * k[c];
  }

  // k4 = f(x + h, y + h*k3)
  computeDerivativeSU3(x + h, xini, tmp, H0, b1_proj, profile,
                       HI_constants, is_antinu, ie, rho, ne, numneu,
                       idata, iglashow, NT_type, nc_factors, k);
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
void __launch_bounds__(128)
evolveKernelImpl(const PhysicsParams params,
                 const PathDeviceData* __restrict__ paths,
                 const double* __restrict__ H0_array,
                 const double* __restrict__ b1_proj,
                 const InteractionDataGPU* __restrict__ interaction_data_ptr,
                 double* __restrict__ d_workspace_corrected,
                 double* __restrict__ d_workspace_sf,
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

  // Per-path slices of the persistent RK4 staging workspace.
  // Layout: [path_idx * nrhos * ne * SU + (rho * ne + ie) * SU + c].
  // This replaces the former per-thread corrected_buf[MAX_PAIRS*SU] +
  // sf_buf[MAX_PAIRS*SU] locals, which were the dominant source of stack
  // spill and register pressure.
  const int path_state_size = nrhos * ne * SU;
  double* path_corrected = d_workspace_corrected + (size_t)path_idx * path_state_size;
  double* path_sf        = d_workspace_sf        + (size_t)path_idx * path_state_size;

  // Shared memory layout for interactions (when enabled):
  // [0 .. nrhos*ne*SU)           : state buffer for cascade computation
  // [nrhos*ne*SU .. + nrhos*3*ne): nc_factors output
  // Dynamic shared memory is allocated in the kernel launch.
  extern __shared__ double smem[];
  const bool do_interactions = params.iinteraction && interaction_data.n_targets > 0;
  double* s_state = smem;                          // [nrhos][ne][SU]
  double* s_nc_factors = smem + nrhos * ne * SU;   // [nrhos][3][ne]

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
    if (do_interactions) {
      // Load current state to shared memory and refresh cascade factors at
      // time x. This is the source-state for sf's k1 substage (and for
      // sh's k1 after the explicit reset between sf and sh below).
      for (int idx = threadIdx.x; idx < nrhos * ne * SU; idx += blockDim.x)
        s_state[idx] = my_state[idx];
      __syncthreads();

      refreshCascadeFactorsSU3(ne, nrhos, NFLV, x, xini,
                                H0_array, b1_proj, interaction_data,
                                path.profile, params.tauregeneration,
                                params.iglashow, params.NT_type,
                                s_state, s_nc_factors);
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
      // ---- Interaction path: substage-refresh RK4 (Hypothesis 2 / Option B).
      //
      //      The CPU's adaptive RKF45 stepper calls PreDerive at *every* RK
      //      substage, so nc_factors / tau_lep_decays / tau_hadlep_decays /
      //      gr_factors are recomputed from the substage state each time
      //      the RHS is evaluated. Option F (commit 1f619ec) restored
      //      Richardson consistency by refreshing once at x+h/2, but it
      //      still froze the cascade source factors across all 4 RK4 sub-
      //      stages of every step. For a two-stage cascade (CC produces
      //      tau, decay produces secondaries) that frozen-source bias is
      //      the dominant remaining error on Test 2 (tau regen, ~10% rel,
      //      abs error ~0.40 unaffected by Option F).
      //
      //      Here we march sf (full step h) and sh (two half-steps h/2)
      //      substage-by-substage, refreshing s_nc_factors cooperatively
      //      from the substage state at each substage time. Both branches
      //      use the same algorithm, so Richardson consistency is preserved.
      //      Cost: ~6x the cascade refreshes per macro step relative to
      //      Option F's 2x; correctness over perf, Perf #2/#3 will reclaim
      //      throughput.

      // Per-thread storage for the at-most MAX_INT_PAIRS (ie, rho) slots
      // each thread owns. With EVOLVE_THREADS=128 and nrhos<=2, MAX_INT_PAIRS=4
      // covers ne up to 256 (well above the ne<=40 used by interaction tests).
      // y_init/acc/sf_local/st_local hold per-slot working state across the
      // cooperative substage syncs.
      constexpr int MAX_INT_PAIRS = 4;

      // Build the slot list for this thread once. n_slots is the number of
      // (ie, rho) pairs actually owned.
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
      // Compile-time/runtime safety: if a thread tried to own more slots
      // than MAX_INT_PAIRS we'd silently drop work. Trip a printf so it
      // shows up in test output rather than producing wrong physics.
      if (n_slots > MAX_INT_PAIRS) {
        if (threadIdx.x == 0 && blockIdx.x == 0 && step_count == 0) {
          printf("ERROR: substage-refresh kernel needs MAX_INT_PAIRS>=%d "
                 "(ne=%d, nrhos=%d, blockDim=%d). Falling back will give "
                 "wrong results.\n", n_slots, ne, nrhos, blockDim.x);
        }
        n_slots = MAX_INT_PAIRS;
      }

      // Lambda-free helpers — these are __device__ functions with explicit
      // captures, since CUDA device lambdas with complex captures have
      // historically silently produced wrong results in this kernel.

      // RK4 substage coefficients — classic 4-stage tableau:
      //   k1 = f(x,         y),                  acc += k1/6,   tmp = y + h/2 * k1
      //   k2 = f(x + h/2,   y + h/2 * k1),       acc += k2/3,   tmp = y + h/2 * k2
      //   k3 = f(x + h/2,   y + h/2 * k2),       acc += k3/3,   tmp = y + h * k3
      //   k4 = f(x + h,     y + h * k3),         acc += k4/6,   y_out = y + h * acc

      // Stage substep sub-times (relative to substep start) for cascade refresh:
      //   substage 0 (k1) is evaluated at substep start.
      //   substages 1,2 (k2,k3) are evaluated at substep start + h_sub/2.
      //   substage 3 (k4) is evaluated at substep start + h_sub.

      // Save the start-of-step nc_factors values before sf consumes them.
      // sh's first substage (k1) evaluates at the start of the macro step
      // and must use the *same* initial cascade source as sf's k1, otherwise
      // sf and sh would solve different ODEs and Richardson would be biased.
      //
      // Easiest implementation: when starting sh we recompute s_nc_factors
      // from my_state at time x. That cost is amortized into the existing
      // sh substage refreshes; we just have to be sure to do it.

      // ===== sf path: one full RK4 step of size h_try =====
      // Prime s_state with current my_state (it already holds my_state values
      // from the start-of-step refresh above). sf substep starts at time x.
      //
      // For each substage we:
      //   1. Each thread computes k for its slots using s_nc_factors (refreshed
      //      at the substage start time by the previous iteration / the start-
      //      of-step block).
      //   2. Each thread updates its accumulator and computes substage tmp,
      //      writing tmp into s_state[slot] for the next refresh.
      //   3. __syncthreads() to publish s_state to all threads.
      //   4. Cooperatively refresh s_nc_factors at the next substage's time.
      //   5. __syncthreads() to publish refreshed s_nc_factors.
      //
      // After substage 4 (k4), each thread combines acc into y_out for sf
      // and writes to path_sf.

      double y_init_sf[MAX_INT_PAIRS][SU];
      double acc_sf[MAX_INT_PAIRS][SU];

      for (int s = 0; s < n_slots; s++) {
        const double* state_ptr = my_state + (slot_rho[s] * ne + slot_ie[s]) * SU;
        #pragma unroll
        for (int c = 0; c < SU; c++) y_init_sf[s][c] = state_ptr[c];
        #pragma unroll
        for (int c = 0; c < SU; c++) acc_sf[s][c] = 0.0;
      }

      // RK4 substage marching for sf.
      // weights[k] is the Butcher RK4 coefficient on k_k for the y_out sum.
      // a[k]      is the Butcher RK4 coefficient on k_k for substage k+1's tmp.
      const double rk_weight[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
      const double rk_a[4]      = {0.5,     0.5,     1.0,     0.0    };  // a[3] unused
      const double rk_t[4]      = {0.0,     0.5,     0.5,     1.0    };  // substage time as fraction of h

      // y_next holds the next substage's evaluation state per slot
      // (the "tmp = y + a*h*k" produced by the previous substage). It also
      // gets staged to s_state[slot] before the cooperative refresh so other
      // threads' refreshes see it. Keeping a per-thread copy avoids any
      // shared-memory read-after-write hazard between writing to s_state
      // and reading it back inside this same thread.
      double y_next_sf[MAX_INT_PAIRS][SU];

      // sf substages (4 stages over size h_try).
      // s_nc_factors entering substage 0 was refreshed from my_state at x.
      // We refresh between substages so the next substage sees the right factors.
      for (int sub = 0; sub < 4; sub++) {
        double x_sub = x + rk_t[sub] * h_try;

        // Each thread: compute k for its slots from current s_nc_factors,
        // update acc, compute next substage's tmp into s_state and y_next.
        for (int s = 0; s < n_slots; s++) {
          int ie = slot_ie[s];
          int rho = slot_rho[s];
          const double* H0 = H0_array + (rho * ne + ie) * SU;
          const double* proj = b1_proj + rho * NFLV * SU;
          bool is_antinu = ((rho == 1) && params.NT_type == 3)
                         || (params.NT_type == 2);

          double y_eval[SU];
          if (sub == 0) {
            #pragma unroll
            for (int c = 0; c < SU; c++) y_eval[c] = y_init_sf[s][c];
          } else {
            #pragma unroll
            for (int c = 0; c < SU; c++) y_eval[c] = y_next_sf[s][c];
          }

          double k[SU];
          computeDerivativeSU3(x_sub, xini, y_eval, H0, proj, path.profile,
                               params.HI_constants, is_antinu,
                               ie, rho, ne, NFLV,
                               interaction_data, params.iglashow, params.NT_type,
                               s_nc_factors, k);

          #pragma unroll
          for (int c = 0; c < SU; c++)
            acc_sf[s][c] += rk_weight[sub] * k[c];

          if (sub < 3) {
            double a = rk_a[sub] * h_try;
            double* s_ptr = s_state + (rho * ne + ie) * SU;
            #pragma unroll
            for (int c = 0; c < SU; c++) {
              double tmp = y_init_sf[s][c] + a * k[c];
              y_next_sf[s][c] = tmp;
              s_ptr[c] = tmp;
            }
          }
        }

        if (sub < 3) {
          __syncthreads();
          // Refresh cascade factors at the next substage's time.
          double x_next = x + rk_t[sub + 1] * h_try;
          refreshCascadeFactorsSU3(ne, nrhos, NFLV, x_next, xini,
                                    H0_array, b1_proj, interaction_data,
                                    path.profile, params.tauregeneration,
                                    params.iglashow, params.NT_type,
                                    s_state, s_nc_factors);
        }
      }

      // Write sf result to path_sf workspace.
      for (int s = 0; s < n_slots; s++) {
        int ie = slot_ie[s];
        int rho = slot_rho[s];
        double* sf_ptr = path_sf + (rho * ne + ie) * SU;
        #pragma unroll
        for (int c = 0; c < SU; c++)
          sf_ptr[c] = y_init_sf[s][c] + h_try * acc_sf[s][c];
      }
      __syncthreads();

      // ===== sh path: two RK4 sub-half-steps of size h_try/2 each =====
      // Reset s_state to my_state and recompute s_nc_factors at time x so
      // sh's k1 sees the same initial cascade source as sf's k1 did.
      for (int idx = threadIdx.x; idx < nrhos * ne * SU; idx += blockDim.x)
        s_state[idx] = my_state[idx];
      __syncthreads();

      refreshCascadeFactorsSU3(ne, nrhos, NFLV, x, xini,
                                H0_array, b1_proj, interaction_data,
                                path.profile, params.tauregeneration,
                                params.iglashow, params.NT_type,
                                s_state, s_nc_factors);

      double y_init_sh[MAX_INT_PAIRS][SU];
      double acc_sh[MAX_INT_PAIRS][SU];

      // First sub-half-step: starts from my_state at time x, advances to x + h/2.
      double h_half = h_try * 0.5;

      for (int s = 0; s < n_slots; s++) {
        const double* state_ptr = my_state + (slot_rho[s] * ne + slot_ie[s]) * SU;
        #pragma unroll
        for (int c = 0; c < SU; c++) y_init_sh[s][c] = state_ptr[c];
        #pragma unroll
        for (int c = 0; c < SU; c++) acc_sh[s][c] = 0.0;
      }

      double y_next_sh[MAX_INT_PAIRS][SU];

      for (int sub = 0; sub < 4; sub++) {
        double x_sub = x + rk_t[sub] * h_half;

        for (int s = 0; s < n_slots; s++) {
          int ie = slot_ie[s];
          int rho = slot_rho[s];
          const double* H0 = H0_array + (rho * ne + ie) * SU;
          const double* proj = b1_proj + rho * NFLV * SU;
          bool is_antinu = ((rho == 1) && params.NT_type == 3)
                         || (params.NT_type == 2);

          double y_eval[SU];
          if (sub == 0) {
            #pragma unroll
            for (int c = 0; c < SU; c++) y_eval[c] = y_init_sh[s][c];
          } else {
            #pragma unroll
            for (int c = 0; c < SU; c++) y_eval[c] = y_next_sh[s][c];
          }

          double k[SU];
          computeDerivativeSU3(x_sub, xini, y_eval, H0, proj, path.profile,
                               params.HI_constants, is_antinu,
                               ie, rho, ne, NFLV,
                               interaction_data, params.iglashow, params.NT_type,
                               s_nc_factors, k);

          #pragma unroll
          for (int c = 0; c < SU; c++)
            acc_sh[s][c] += rk_weight[sub] * k[c];

          if (sub < 3) {
            double a = rk_a[sub] * h_half;
            double* s_ptr = s_state + (rho * ne + ie) * SU;
            #pragma unroll
            for (int c = 0; c < SU; c++) {
              double tmp = y_init_sh[s][c] + a * k[c];
              y_next_sh[s][c] = tmp;
              s_ptr[c] = tmp;
            }
          }
        }

        if (sub < 3) {
          __syncthreads();
          double x_next = x + rk_t[sub + 1] * h_half;
          refreshCascadeFactorsSU3(ne, nrhos, NFLV, x_next, xini,
                                    H0_array, b1_proj, interaction_data,
                                    path.profile, params.tauregeneration,
                                    params.iglashow, params.NT_type,
                                    s_state, s_nc_factors);
        }
      }

      // st = result of first sub-half-step (kept per-thread; also published
      // to s_state for the second sub-half-step's k1 cascade refresh).
      double st_local[MAX_INT_PAIRS][SU];
      for (int s = 0; s < n_slots; s++) {
        int ie = slot_ie[s];
        int rho = slot_rho[s];
        #pragma unroll
        for (int c = 0; c < SU; c++)
          st_local[s][c] = y_init_sh[s][c] + h_half * acc_sh[s][c];
        double* s_ptr = s_state + (rho * ne + ie) * SU;
        #pragma unroll
        for (int c = 0; c < SU; c++) s_ptr[c] = st_local[s][c];
      }
      __syncthreads();

      // Refresh cascade factors at time x + h_half using the half-step result.
      refreshCascadeFactorsSU3(ne, nrhos, NFLV, x + h_half, xini,
                                H0_array, b1_proj, interaction_data,
                                path.profile, params.tauregeneration,
                                params.iglashow, params.NT_type,
                                s_state, s_nc_factors);

      // Second sub-half-step: starts from st at time x + h_half, advances to x + h.
      // Reuse y_init_sh / acc_sh as the working storage for the second half.
      for (int s = 0; s < n_slots; s++) {
        #pragma unroll
        for (int c = 0; c < SU; c++) y_init_sh[s][c] = st_local[s][c];
        #pragma unroll
        for (int c = 0; c < SU; c++) acc_sh[s][c] = 0.0;
      }

      for (int sub = 0; sub < 4; sub++) {
        double x_sub = (x + h_half) + rk_t[sub] * h_half;

        for (int s = 0; s < n_slots; s++) {
          int ie = slot_ie[s];
          int rho = slot_rho[s];
          const double* H0 = H0_array + (rho * ne + ie) * SU;
          const double* proj = b1_proj + rho * NFLV * SU;
          bool is_antinu = ((rho == 1) && params.NT_type == 3)
                         || (params.NT_type == 2);

          double y_eval[SU];
          if (sub == 0) {
            #pragma unroll
            for (int c = 0; c < SU; c++) y_eval[c] = y_init_sh[s][c];
          } else {
            #pragma unroll
            for (int c = 0; c < SU; c++) y_eval[c] = y_next_sh[s][c];
          }

          double k[SU];
          computeDerivativeSU3(x_sub, xini, y_eval, H0, proj, path.profile,
                               params.HI_constants, is_antinu,
                               ie, rho, ne, NFLV,
                               interaction_data, params.iglashow, params.NT_type,
                               s_nc_factors, k);

          #pragma unroll
          for (int c = 0; c < SU; c++)
            acc_sh[s][c] += rk_weight[sub] * k[c];

          if (sub < 3) {
            double a = rk_a[sub] * h_half;
            double* s_ptr = s_state + (rho * ne + ie) * SU;
            #pragma unroll
            for (int c = 0; c < SU; c++) {
              double tmp = y_init_sh[s][c] + a * k[c];
              y_next_sh[s][c] = tmp;
              s_ptr[c] = tmp;
            }
          }
        }

        if (sub < 3) {
          __syncthreads();
          double x_next = (x + h_half) + rk_t[sub + 1] * h_half;
          refreshCascadeFactorsSU3(ne, nrhos, NFLV, x_next, xini,
                                    H0_array, b1_proj, interaction_data,
                                    path.profile, params.tauregeneration,
                                    params.iglashow, params.NT_type,
                                    s_state, s_nc_factors);
        }
      }

      // Combine into corrected output via Richardson extrapolation.
      // sh = y_init_sh + h_half * acc_sh, and sf is in path_sf.
      for (int s = 0; s < n_slots; s++) {
        int ie = slot_ie[s];
        int rho = slot_rho[s];
        double sh[SU];
        #pragma unroll
        for (int c = 0; c < SU; c++)
          sh[c] = y_init_sh[s][c] + h_half * acc_sh[s][c];

        const double* sf_ptr = path_sf + (rho * ne + ie) * SU;
        const double* y0 = my_state + (rho * ne + ie) * SU;
        double* corr = path_corrected + (rho * ne + ie) * SU;

        #pragma unroll
        for (int c = 0; c < SU; c++) {
          corr[c] = sh[c] + (sh[c] - sf_ptr[c]) / 15.0;
          double err = fabs(sf_ptr[c] - sh[c]) / 15.0;
          double scale = solver_config.abs_error +
                         solver_config.rel_error * fmax(fabs(y0[c]), fabs(sh[c]));
          if (scale > 0.0)
            local_max_err = fmax(local_max_err, err / scale);
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
      // Phase 3: Write corrected states to global memory
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
                  const SolverConfig& solver_config,
                  double* d_states,
                  int n_paths, int numneu,
                  cudaStream_t stream) {
  if (n_paths == 0) return;

  constexpr int threads = EVOLVE_THREADS;

  // Compute shared memory for interaction cascade computation
  // Layout: s_state[nrhos*ne*SU] + s_nc_factors[nrhos*3*ne]
  size_t shared_bytes = 0;
  if (params.iinteraction && d_interaction_data && d_interaction_data->n_targets > 0) {
    int su_size = numneu * numneu;
    shared_bytes = (size_t)(params.nrhos * params.ne * su_size   // state
                          + params.nrhos * 3 * params.ne)        // nc_factors
                   * sizeof(double);
  }

  // Upload InteractionDataGPU struct to device if needed
  // (the struct itself contains device pointers, but it needs to be device-resident
  //  for the kernel to read it via a device pointer)
  InteractionDataGPU* d_idata_on_device = nullptr;
  if (d_interaction_data && d_interaction_data->n_targets > 0) {
    NUSQUIDS_CUDA_CHECK(cudaMalloc(&d_idata_on_device, sizeof(InteractionDataGPU)));
    NUSQUIDS_CUDA_CHECK(cudaMemcpyAsync(d_idata_on_device, d_interaction_data,
                                         sizeof(InteractionDataGPU),
                                         cudaMemcpyHostToDevice, stream));
  }

  // Clear any stale CUDA errors (e.g., from cudaDeviceSetLimit on MIG)
  cudaGetLastError();

  switch (numneu) {
    case 3:
      evolveKernelImpl<3><<<n_paths, threads, shared_bytes, stream>>>(
        params, d_paths, d_H0_array, d_b1_proj, d_idata_on_device,
        d_workspace_corrected, d_workspace_sf,
        solver_config, d_states, n_paths);
      break;
    case 4:
      evolveKernelImpl<4><<<n_paths, threads, shared_bytes, stream>>>(
        params, d_paths, d_H0_array, d_b1_proj, d_idata_on_device,
        d_workspace_corrected, d_workspace_sf,
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
