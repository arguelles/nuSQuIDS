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
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Forward declarations — no CUDA includes in this header
namespace nusquids {

class nuSQUIDS;
template<typename,typename> class nuSQUIDSAtm;
class EarthAtm;
class Body;

/// \brief GPU kernel launch-configuration constants — mirrored in
/// `include/nuSQuIDS/cuda/detail/kernels.cuh`. Duplicated here (pure-host
/// header, no CUDA includes) so `nuSQUIDS::EvolveState` can enforce the
/// same host-side preconditions without pulling in nvcc-only headers. The
/// kernels.cuh copies are the source of truth for the device kernel;
/// values must match — a static_assert in kernels.cu ties them together.
static constexpr int GPU_EVOLVE_THREADS = 128;
static constexpr int GPU_MAX_INT_PAIRS = 4;

/// \brief Return true iff the GPU pure-oscillation integrator will lock up
/// on the given configuration and the CPU integrator should be used instead.
///
/// The GPU pure-osc code path uses RK4 + Richardson step doubling, which
/// cannot resolve fast oscillations superposed on a matter potential — the
/// stepper collapses h to h_min and burns through max_steps before finishing.
/// This is empirically what happened at Emin=10 MeV through Earth matter, and
/// at Emin=100 MeV through ~500 km of rock. Vacuum works because there is no
/// matter Hamiltonian to interfere with the RK4 truncation error. The
/// interacting substage-refresh path has independent step control and handles
/// low energies fine. So the rule is:
///
///   fallback ⇔ (iinteraction == false) AND (matter density > 0 somewhere)
///
/// Higher energies (Emin ≳ 10 GeV) with matter + osc-only actually work fine
/// too, but a conservative rule avoids surprising failures for users. This
/// keeps atmospheric-oscillation grids (the paper's 200×200 workload) on
/// GPU because they always enable interactions.
inline bool GPUShouldFallBackToCPU(bool iinteraction, bool matter_present) {
  return !iinteraction && matter_present;
}

/// \brief Throw std::runtime_error with a clear message if the given
/// (numneu, ne, nrhos, iinteraction, n_targets) configuration is outside
/// what the GPU backend supports.
///
/// Called by both `nuSQUIDS::EvolveState` and `nuSQUIDSAtm::EvolveState`
/// before dispatching to `CUDABackend::Evolve`. Identical checks live in
/// `launchEvolve` (kernels.cu) as defence in depth, but experience with
/// nvcc-12.4.1 + `--use_fast_math` + `-O3` showed those throws could be
/// optimised away silently in the kernel-launcher TU. Enforcing here
/// (a pure-host TU compiled by g++/clang++) guarantees the check runs.
inline void CheckGPUPreconditions(int numneu, int ne, int nrhos,
                                  bool iinteraction, int n_targets) {
  if (numneu != 3) {
    std::ostringstream msg;
    msg << "nuSQUIDS::EvolveState: SU(" << numneu
        << ") evolution is not implemented on the GPU backend. "
        << "Only numneu=3 is supported; please run on the CPU backend "
        << "(Set_Backend(Backend::cpu)) or set numneu=3.";
    throw std::runtime_error(msg.str());
  }

  if (iinteraction && n_targets > 0) {
    const int int_pairs_per_thread =
        nrhos * ((ne + GPU_EVOLVE_THREADS - 1) / GPU_EVOLVE_THREADS);
    if (int_pairs_per_thread > GPU_MAX_INT_PAIRS) {
      std::ostringstream msg;
      msg << "nuSQUIDS::EvolveState: GPU substage-refresh requires "
          << "MAX_INT_PAIRS>=" << int_pairs_per_thread
          << " (ne=" << ne
          << ", nrhos=" << nrhos
          << ", blockDim=" << GPU_EVOLVE_THREADS
          << ", current MAX_INT_PAIRS=" << GPU_MAX_INT_PAIRS << ").";
      throw std::runtime_error(msg.str());
    }
  }
}

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
  int rounded_ne;              ///< ne rounded up to preferred_alignment (stride for last dim)

  bool has_glashow;       ///< Whether Glashow resonance data is present
  bool has_tau_regen;     ///< Whether tau regeneration data is present

  InteractionDataHost() :
    n_targets(0), nrhos(0), numneu(0), ne(0), rounded_ne(0),
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

  // Electron number density along the track (natural units, eV³),
  // precomputed on CPU using the same logic as nuSQUIDS::UpdateInteractions
  // (handles isoscalar, p/n, body-composition, and nuclear-XS branches in
  // squids::Const natural units). Used by Glashow on the GPU. Empty when
  // Glashow is disabled.
  std::vector<double> num_e_vals;                    ///< [n_density_samples]
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
