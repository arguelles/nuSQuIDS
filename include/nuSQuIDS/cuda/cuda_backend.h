/******************************************************************************
*    This program is free software: you can redistribute it and/or modify     *
*   it under the terms of the GNU General Public License as published by      *
*   the Free Software Foundation, either version 3 of the License, or         *
*   (at your option) any later version.                                       *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License         *
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
******************************************************************************/

#ifndef NUSQUIDS_CUDA_BACKEND_H
#define NUSQUIDS_CUDA_BACKEND_H

#include <memory>
#include <vector>
#include <string>

// Forward declarations — no CUDA includes in this header
namespace nusquids {

class nuSQUIDS;
template<typename,typename> class nuSQUIDSAtm;
class EarthAtm;
class Body;

/// \brief Host-side interaction data for GPU upload.
///
/// Holds flat pointers to the cross-section tables from InteractionStructure.
/// Filled by the caller (nuSQUIDS) and passed to CUDABackend::Evolve().
/// No CUDA types — safe to include from any translation unit.
struct InteractionDataHost {
  int n_targets;     ///< Number of nuclear target types
  int nrhos;         ///< Number of density matrix equations
  int numneu;        ///< Number of neutrino flavors
  int ne;            ///< Number of energy nodes

  const double* sigma_CC;      ///< Total CC cross section [n_targets * nrhos * numneu * ne]
  const double* sigma_NC;      ///< Total NC cross section [n_targets * nrhos * numneu * ne]
  const double* dNdE_CC;       ///< Differential CC xsec  [n_targets * nrhos * numneu * ne * ne]
  const double* dNdE_NC;       ///< Differential NC xsec  [n_targets * nrhos * numneu * ne * ne]
  const double* sigma_GR;      ///< Glashow cross section [ne] (may be nullptr)
  const double* dNdE_GR;       ///< Glashow differential  [ne * ne] (may be nullptr)
  const double* dNdE_tau_all;  ///< Tau decay → all       [nrhos * ne * ne] (may be nullptr)
  const double* dNdE_tau_lep;  ///< Tau decay → leptons    [nrhos * ne * ne] (may be nullptr)
  const double* energies;      ///< Energy nodes [ne]
  const double* delE;          ///< Energy bin widths [ne-1]

  bool has_glashow;       ///< Whether Glashow resonance data is present
  bool has_tau_regen;     ///< Whether tau regeneration data is present

  InteractionDataHost() :
    n_targets(0), nrhos(0), numneu(0), ne(0),
    sigma_CC(nullptr), sigma_NC(nullptr),
    dNdE_CC(nullptr), dNdE_NC(nullptr),
    sigma_GR(nullptr), dNdE_GR(nullptr),
    dNdE_tau_all(nullptr), dNdE_tau_lep(nullptr),
    energies(nullptr), delE(nullptr),
    has_glashow(false), has_tau_regen(false) {}
};

/// \brief Parameters needed for GPU propagation of a single path.
struct GPUPathData {
  double xini;         ///< Start position along track
  double xend;         ///< End position along track
  double time_offset;  ///< Time offset for evolution
  int n_density_samples;  ///< Number of density profile samples
  std::vector<double> density_x;    ///< Sample positions
  std::vector<double> density_vals; ///< Density at each sample (g/cm³ for HI computation)
  std::vector<double> ye_vals;      ///< Electron fraction at each sample

  // Target number densities along the track (natural units, eV³)
  // Precomputed on CPU using GetTargetNumberDensities() from squids::Const.
  // n_targets arrays, each of size n_density_samples.
  int n_targets;                                     ///< Number of target types
  std::vector<std::vector<double>> target_ndens;     ///< [n_targets][n_density_samples]
};

/// \brief CUDA GPU backend for nuSQuIDS propagation.
///
/// This class provides GPU-accelerated neutrino propagation using CUDA.
/// It uses the pimpl pattern so that no CUDA headers are needed by users.
/// The GPU backend works with raw state arrays to avoid needing friend access
/// to nuSQUIDS internals.
class CUDABackend {
public:
  struct Config {
    std::vector<int> device_ids;  ///< GPU device IDs to use (empty = auto-detect all)
    int batch_size_limit;         ///< Max paths per GPU batch
    double rel_error;             ///< Relative error tolerance for ODE solver
    double abs_error;             ///< Absolute error tolerance for ODE solver
    Config() : batch_size_limit(2000), rel_error(1e-6), abs_error(1e-6) {}
  };

  CUDABackend(Config config = Config());
  ~CUDABackend();

  // Non-copyable, movable
  CUDABackend(const CUDABackend&) = delete;
  CUDABackend& operator=(const CUDABackend&) = delete;
  CUDABackend(CUDABackend&&) noexcept;
  CUDABackend& operator=(CUDABackend&&) noexcept;

  /// \brief Evolve all paths on GPU.
  ///
  /// \param states     Flat array of SU_vector components for all paths.
  ///                   Layout: [n_paths][nrhos][ne][su_size].
  ///                   Modified in-place with evolved states.
  /// \param paths      Per-path data (body profile, track endpoints).
  /// \param H0_array   Vacuum Hamiltonian components.
  ///                   Layout: [nrhos][ne][su_size].
  /// \param b1_proj    Flavor projectors.
  ///                   Layout: [nrhos][numneu][su_size].
  /// \param n_paths    Number of paths (zenith angles).
  /// \param ne         Number of energy nodes.
  /// \param nrhos      Number of density matrix types (1 or 2).
  /// \param numneu     Number of neutrino flavors.
  /// \param HI_constants  Matter potential constant: sqrt(2)*GF*Na/cm^3 [eV].
  /// \param NT_type    Neutrino type: 1=neutrino, 2=antineutrino, 3=both.
  void Evolve(double* states,
              const std::vector<GPUPathData>& paths,
              const double* H0_array,
              const double* b1_proj,
              int n_paths, int ne, int nrhos, int numneu,
              double HI_constants, int NT_type,
              double rel_error = 0.0, double abs_error = 0.0,
              const InteractionDataHost* interaction_data = nullptr);

  /// \brief Runtime GPU detection — returns true if any CUDA GPU is available.
  static bool IsAvailable();

  /// \brief Returns the number of available CUDA GPUs.
  static int DeviceCount();

  /// \brief Returns a string describing the available CUDA device(s).
  static std::string DeviceInfo();

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace nusquids

#endif // NUSQUIDS_CUDA_BACKEND_H
