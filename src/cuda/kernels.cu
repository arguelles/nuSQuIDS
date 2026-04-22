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
// Mirrors the CPU logic in src/nuSQuIDS.cpp UpdateInteractions(),
// without the composition/nuclear-XS branch (not yet plumbed to GPU):
//   - n_targets == 1 (isoscalar nucleon): n_e = N_nuc * ye
//   - n_targets >= 2 (proton/neutron):    n_e = N_p  (neutral matter)
// ============================================================

__device__ __forceinline__
double electronNumberDensitySU3(const double* __restrict__ target_ndens,
                                int n_targets, double ye)
{
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
      double num_e = electronNumberDensitySU3(target_ndens, idata.n_targets, ye);
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

      // Decay contributions: tau → neutrinos at energy e1
      if (tau_flux > 0.0 || taubar_flux > 0.0) {
        // dNdE_tau_all[rho][et][e1], dNdE_tau_lep[rho][et][e1] with rounded_ne stride
        int rne = idata.rounded_ne;
        size_t tau_idx_nu    = (0 * ne + et) * rne + e1;
        size_t tau_idx_nubar = (1 * ne + et) * rne + e1;

        tau_hadlep_nu    += tau_flux    * idata.d_dNdE_tau_all[tau_idx_nu];
        tau_lep_nubar    += tau_flux    * idata.d_dNdE_tau_lep[tau_idx_nubar];
        tau_hadlep_nubar += taubar_flux * idata.d_dNdE_tau_all[tau_idx_nubar];
        tau_lep_nu       += taubar_flux * idata.d_dNdE_tau_lep[tau_idx_nu];
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

  // Compute electron number density from the same target splines used by
  // NC/CC absorption. In neutral matter n_e = n_p (proton target) when
  // n_targets>=2; for isoscalar nucleon targets, n_e = N_nuc * ye.
  double target_ndens[MAX_TARGETS];
  for (int t = 0; t < idata.n_targets && t < MAX_TARGETS; t++)
    target_ndens[t] = evaluateTargetFraction(profile, t, x_eval);
  double ye = evaluateYe(profile, x_eval);
  double num_e = electronNumberDensitySU3(target_ndens, idata.n_targets, ye);

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
                   iglashow, NT_type, ye_x, invlen);

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

  // Buffer for corrected states: each thread handles at most
  // ceil(ne/blockDim) * nrhos pairs, each with SU components.
  // The invariant nrhos * ceil(ne / EVOLVE_THREADS) <= MAX_PAIRS is
  // checked on the host side in launchEvolve() before kernel launch.
  double corrected_buf[MAX_PAIRS * SU];

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
      // Load current state to shared memory
      for (int idx = threadIdx.x; idx < nrhos * ne * SU; idx += blockDim.x)
        s_state[idx] = my_state[idx];
      __syncthreads();

      // Compute NC cascade factors cooperatively
      // Target number densities are interpolated from profile splines
      // (precomputed on CPU in natural units using squids::Const)
      computeNCCascadeSU3(ne, nrhos, NFLV, x, xini,
                          H0_array, b1_proj, interaction_data, path.profile,
                          s_state, s_nc_factors);
      __syncthreads();


      // Tau regeneration (adds to nc_factors in-place).
      // Uses profile-based target number densities via evaluateTargetFraction.
      if (params.tauregeneration && interaction_data.d_dNdE_tau_all) {
        computeTauRegenSU3(ne, nrhos, NFLV, x, xini,
                           H0_array, b1_proj, interaction_data,
                           path.profile, s_state, s_nc_factors);
        __syncthreads();
      }

      // Glashow resonance (adds to nc_factors in-place).
      // Uses profile-based electron number density via evaluateTargetFraction.
      if (params.iglashow && interaction_data.d_sigma_GR) {
        computeGlashowCascadeSU3(ne, nrhos, NFLV, x, xini,
                                  H0_array, b1_proj, interaction_data,
                                  path.profile, params.NT_type,
                                  s_state, s_nc_factors);
        __syncthreads();
      }
    }

    double local_max_err = 0.0;
    int pair_idx = 0;

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
    } else {
      // ---- Interaction path: refresh nc_factors between half-steps ----
      // Buffer for full-step results (needed for Richardson at end)
      double sf_buf[MAX_PAIRS * SU];

      // Pass 1: Full step + first half-step using initial nc_factors.
      // Write first half-step results (st) to shared memory for cascade refresh.
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

          // Full step → sf_buf
          double* sf = sf_buf + pair_idx * SU;
          rk4StepSU3_interacting(y, x, h_try, xini, H0, proj, path.profile,
                    params.HI_constants, is_antinu, ie, rho, ne, NFLV,
                    interaction_data, params.iglashow, params.NT_type,
                    s_nc_factors, sf);

          // First half-step → s_state (shared memory) for cascade refresh
          double st[SU];
          rk4StepSU3_interacting(y, x, h_try * 0.5, xini, H0, proj, path.profile,
                    params.HI_constants, is_antinu, ie, rho, ne, NFLV,
                    interaction_data, params.iglashow, params.NT_type,
                    s_nc_factors, st);

          // Store st to shared memory so cascade can read all energies
          double* s_st = s_state + (rho * ne + ie) * SU;
          #pragma unroll
          for (int c = 0; c < SU; c++) s_st[c] = st[c];

          // Also store st to corrected_buf temporarily (reused in Pass 2)
          double* st_buf = corrected_buf + pair_idx * SU;
          #pragma unroll
          for (int c = 0; c < SU; c++) st_buf[c] = st[c];

          pair_idx++;
        }
      }
      __syncthreads();

      // Refresh nc_factors from mid-step state at position x + h_try/2
      double x_mid = x + h_try * 0.5;
      computeNCCascadeSU3(ne, nrhos, NFLV, x_mid, xini,
                          H0_array, b1_proj, interaction_data, path.profile,
                          s_state, s_nc_factors);
      __syncthreads();

      if (params.tauregeneration && interaction_data.d_dNdE_tau_all) {
        computeTauRegenSU3(ne, nrhos, NFLV, x_mid, xini,
                           H0_array, b1_proj, interaction_data,
                           path.profile, s_state, s_nc_factors);
        __syncthreads();
      }

      if (params.iglashow && interaction_data.d_sigma_GR) {
        computeGlashowCascadeSU3(ne, nrhos, NFLV, x_mid, xini,
                                  H0_array, b1_proj, interaction_data,
                                  path.profile, params.NT_type,
                                  s_state, s_nc_factors);
        __syncthreads();
      }

      // Pass 2: Second half-step using refreshed nc_factors + Richardson
      pair_idx = 0;
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

          // Read st from corrected_buf (stored in Pass 1)
          double* st = corrected_buf + pair_idx * SU;

          // Second half-step using refreshed nc_factors
          double sh[SU];
          rk4StepSU3_interacting(st, x + h_try * 0.5, h_try * 0.5, xini, H0, proj, path.profile,
                    params.HI_constants, is_antinu, ie, rho, ne, NFLV,
                    interaction_data, params.iglashow, params.NT_type,
                    s_nc_factors, sh);

          // Richardson extrapolation and error estimate
          const double* sf = sf_buf + pair_idx * SU;
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
                  const SolverConfig& solver_config,
                  double* d_states,
                  int n_paths, int numneu,
                  cudaStream_t stream) {
  if (n_paths == 0) return;

  constexpr int threads = EVOLVE_THREADS;

  // Guard: the kernel sizes per-thread RK4 correction buffers for at most
  // MAX_PAIRS (rho, ie) pairs. Silently exceeding this overruns thread-local
  // memory and produces non-local numerical corruption.
  {
    const int pairs_per_thread = params.nrhos * ((params.ne + threads - 1) / threads);
    if (pairs_per_thread > MAX_PAIRS) {
      throw std::runtime_error(
        "nuSQuIDS CUDA backend: evolveKernel per-thread pair count " +
        std::to_string(pairs_per_thread) + " exceeds MAX_PAIRS=" +
        std::to_string(MAX_PAIRS) + " (ne=" + std::to_string(params.ne) +
        ", nrhos=" + std::to_string(params.nrhos) +
        ", threads=" + std::to_string(threads) +
        "). Reduce ne or raise MAX_PAIRS/EVOLVE_THREADS in kernels.cuh.");
    }
  }

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
        params, d_paths, d_H0_array, d_b1_proj,
        d_idata_on_device, solver_config, d_states, n_paths);
      break;
    case 4:
      evolveKernelImpl<4><<<n_paths, threads, shared_bytes, stream>>>(
        params, d_paths, d_H0_array, d_b1_proj,
        d_idata_on_device, solver_config, d_states, n_paths);
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
