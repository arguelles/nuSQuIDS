/******************************************************************************
*   Cross-section data on GPU for nuSQuIDS CUDA backend.
*   Uploads InteractionStructure to GPU once (shared across all paths).
*   Uses SoA layout with energy as the contiguous dimension for coalesced access.
******************************************************************************/

#ifndef NUSQUIDS_CUDA_INTERACTIONS_GPU_CUH
#define NUSQUIDS_CUDA_INTERACTIONS_GPU_CUH

#include <cuda_runtime.h>
#include <vector>
#include "memory.cuh"

namespace nusquids { namespace cuda {

/// GPU-side interaction data layout
/// All arrays have energy as the fast (contiguous) dimension for coalesced access.
struct InteractionDataGPU {
  int n_targets;       ///< Number of nuclear target types
  int nrhos;           ///< Number of density matrix equations
  int numneu;          ///< Number of neutrino flavors
  int ne;              ///< Number of energy nodes

  // Cross-section tables on device
  // Layout: [target][rho][flavor][ie_in][ie_out] for dNdE
  //         [target][rho][flavor][ie] for sigma

  double* d_dNdE_CC;   ///< Differential CC cross section: [targets * nrhos * numneu * ne * ne]
  double* d_dNdE_NC;   ///< Differential NC cross section
  double* d_sigma_CC;  ///< Total CC cross section: [targets * nrhos * numneu * ne]
  double* d_sigma_NC;  ///< Total NC cross section

  // Tau decay spectra
  double* d_dNdE_tau_all; ///< Tau decay to all: [ne * ne]
  double* d_dNdE_tau_lep; ///< Tau decay to leptons: [ne * ne]

  // Glashow resonance
  double* d_sigma_GR;    ///< Glashow cross section: [ne]
  double* d_dNdE_GR;     ///< Glashow differential: [ne * ne]

  // Energy grid
  double* d_energies;    ///< Energy nodes [ne]
  double* d_delE;        ///< Energy bin widths [ne-1]

  // Flavor projectors (initial, before evolution)
  double* d_b0_proj;    ///< b0 projectors: [numneu * su_size]
  double* d_b1_proj;    ///< b1 projectors: [nrhos * numneu * su_size]

  /// Total device memory used (bytes)
  size_t total_bytes;
};

/// Upload InteractionStructure data to GPU
/// \param int_struct Pointer to the CPU-side InteractionStructure
/// \param ne Number of energy nodes
/// \param nrhos Number of density matrix equations
/// \param numneu Number of neutrino flavors
/// \param stream CUDA stream for async operations
/// \return Populated InteractionDataGPU with device pointers
InteractionDataGPU uploadInteractionData(const void* int_struct,
                                         int ne, int nrhos, int numneu,
                                         const double* energies,
                                         const double* delE,
                                         const double* b0_proj,
                                         const double* b1_proj,
                                         cudaStream_t stream = 0);

/// Free all device memory in InteractionDataGPU
void freeInteractionData(InteractionDataGPU& data);

/// Helper: compute flat index for dNdE arrays
/// Layout: [target][rho][flavor][ie_in * ne + ie_out]
__device__ __host__ __forceinline__
size_t dNdE_index(int target, int rho, int flavor, int ie_in, int ie_out,
                  int nrhos, int numneu, int ne) {
  return ((((size_t)target * nrhos + rho) * numneu + flavor) * ne + ie_in) * ne + ie_out;
}

/// Helper: compute flat index for sigma arrays
/// Layout: [target][rho][flavor][ie]
__device__ __host__ __forceinline__
size_t sigma_index(int target, int rho, int flavor, int ie,
                   int nrhos, int numneu, int ne) {
  return (((size_t)target * nrhos + rho) * numneu + flavor) * ne + ie;
}

}} // namespace nusquids::cuda

#endif // NUSQUIDS_CUDA_INTERACTIONS_GPU_CUH
