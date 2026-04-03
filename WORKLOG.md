# nuSQuIDS GPU Backend Work Log

## 2025-02-20 — Initial GPU Backend Implementation

**What was done:**
- Implemented CUDA GPU backend for `nuSQUIDSAtm` (atmospheric, multi-path propagation)
- Oscillation-only mode, SU(3) (3 neutrino flavors)
- Adaptive RK4 solver with Richardson extrapolation and PI step-size controller
- Akima spline density profile interpolation on GPU
- Multi-GPU support with round-robin path distribution
- Pimpl pattern to hide CUDA dependencies from public API

**Files created:**
- `include/nuSQuIDS/cuda/cuda_backend.h` — public CUDABackend API
- `include/nuSQuIDS/cuda/detail/body_gpu.cuh` — density profile + Akima splines
- `include/nuSQuIDS/cuda/detail/interactions_gpu.cuh` — cross-section structs (TODO)
- `include/nuSQuIDS/cuda/detail/kernels.cuh` — kernel launch declarations
- `include/nuSQuIDS/cuda/detail/memory.cuh` — RAII GPU memory management
- `include/nuSQuIDS/cuda/detail/physics.cuh` — physics parameter structs
- `include/nuSQuIDS/cuda/detail/propagator.cuh` — per-GPU propagator class
- `include/nuSQuIDS/cuda/detail/solver_gpu.cuh` — ODE solver config
- `include/nuSQuIDS/cuda/detail/sumath.cuh` — SU(N) algebra (SU(3) complete)
- `src/cuda/cuda_backend.cu` — GPU detection, path conversion, Akima fitting
- `src/cuda/kernels.cu` — iCommutator, HI computation, RK4 step, main kernel
- `src/cuda/propagator.cu` — memory management, data upload/download, batching

**Files modified:**
- `include/nuSQuIDS/nuSQuIDS.h` — Backend enum, CUDABackend forward decl, nuSQUIDSAtm GPU path
- `configure` — `--enable-cuda` option, nvcc detection, GPU architecture flags
- `Makefile` — CUDA source compilation rules

**Test results (A100):**
- All tests PASS
- Max absolute error: ~1e-6 (GPU vs CPU)
- Max relative error: ~0.1%
- Speedup: up to 40x for large problem sizes (200 zenith × 200 energy)

**Issues resolved:**
- GSL 2.7→2.8 compatibility fix for cluster
- `deploy_cuda.sh` fixes: SQuIDS path, `CXX=g++`, sbatch path
- `.gitignore` fix for `test/env_vars.sh`

**Known limitations:**
- Only SU(3) — SU(4) through SU(6) stubbed out (kernel returns early)
- No interactions on GPU (NC cascading, tau regen, Glashow)
- Cross-section data not uploaded to GPU (`propagator.cu:84` TODO)
- Only `nuSQUIDSAtm` had GPU support (not single-baseline `nuSQUIDS`)

---

## 2025-02-20 — Test Infrastructure

**What was done:**
- Created comprehensive CUDA test (`test/cuda/cuda_test.cpp`)
  - Test 1: GPU detection (IsAvailable, DeviceCount, DeviceInfo)
  - Test 2: Backend enum and Set/Get API
  - Test 3: GPU vs CPU propagation with detailed per-flavor diagnostics
- Created benchmark suite (`test/cuda/benchmark_cuda.cpp`)
  - Various problem sizes: ncz ∈ {5,10,20,40,100,200}, ne ∈ {20,40,60,100,200}
- Created SLURM test runner (`test/cuda/run_cuda_test.sh`)
- Created deployment script (`deploy_cuda.sh`)
- Added Python (pybind11) and Julia bindings for Backend enum and Set/Get

---

## 2026-02-20 — Single-Baseline GPU Support

**What was done:**
- Added `backend_` member and `cuda_backend_` to nuSQUIDS base class
- Added `Set_Backend()`/`Get_Backend()` methods to nuSQUIDS
- Implemented GPU evolution path in `nuSQUIDS::EvolveState()`
  - Multi-energy only (ne > 1); falls back to CPU for single-energy
  - Packages H0, b1_proj, states, density profile for CUDABackend::Evolve()
  - Supports oscillation-only mode (same as nuSQUIDSAtm GPU path)
- Created single-baseline CUDA test (`test/cuda/cuda_single_baseline_test.cpp`)
  - Vacuum GPU vs CPU comparison
  - Earth matter GPU vs CPU comparison
- Updated Python/Julia bindings to expose Set_Backend/Get_Backend on nuSQUIDS
- Updated `test/cuda/run_cuda_test.sh` with new test
- Created project documentation: CLAUDE.md, skills.md, WORKLOG.md

**Files modified:**
- `include/nuSQuIDS/nuSQuIDS.h` — nuSQUIDS class: backend_, cuda_backend_, Set/Get_Backend
- `src/nuSQuIDS.cpp` — GPU path in EvolveState()
- `resources/python/src_pybind/nuSQUIDSpy.h` — Set/Get_Backend for nuSQUIDS
- `resources/julia/src/nuSQuIDSjl.cpp` — Set/Get_Backend for nuSQUIDS
- `test/cuda/run_cuda_test.sh` — compile + run new test

**Files created:**
- `test/cuda/cuda_single_baseline_test.cpp`
- `CLAUDE.md`, `skills.md`, `WORKLOG.md`

---

## 2026-04-02 — Full GPU Interaction Support (Phases 1-6)

**What was done:**
- Phase 1: Cross-section data upload infrastructure
  - `InteractionDataHost` struct in `cuda_backend.h` (public API, no CUDA deps)
  - `uploadInteractionData()` uploads sigma_CC/NC, dNdE_CC/NC, Glashow, tau spectra
  - Multi-GPU safe (each Propagator uploads its own copy)
- Phase 2: Per-step invlen computation kernel
  - `computeInvlenSU3()` — inverse interaction lengths from cross-sections × density
  - ye-based target fractions for proton/neutron (standard case)
- Phase 3: GammaRho absorption kernel
  - `antiCommutatorSU3()` — translated from SQuIDS AnticonmutatorSU3.txt
  - `evolveProjectorsSU3()` — factored out projector evolution for reuse
  - `computeGammaRhoSU3()` — absorption as weighted evolved projector sum
  - `rk4StepSU3_interacting()` — full interacting RK4: dρ/dt = i[ρ,HI] - {Γ,ρ} + F
- Phase 4: NC cascade kernel (O(ne²))
  - `computeNCCascadeSU3()` — cooperative shared-memory cascade computation
  - Lagged approximation: nc_factors computed once per adaptive step
  - Dynamic shared memory: [nrhos×ne×SU + nrhos×3×ne] doubles
- Phase 5: Tau regeneration + Glashow resonance
  - `computeTauRegenSU3()` — two-pass O(ne²) tau production + decay
  - `computeGlashowCascadeSU3()` — electron antineutrino Glashow cascade
- Phase 6: Full wiring
  - Removed `!iinteraction` guard from GPU path in nuSQUIDS::EvolveState()
  - Added InteractionDataHost extraction from int_struct in both nuSQUIDS and nuSQUIDSAtm
  - Kernel dispatches to interacting vs oscillation-only RK4 based on params.iinteraction

**Files modified:**
- `include/nuSQuIDS/cuda/cuda_backend.h` — InteractionDataHost struct, extended Evolve() API
- `include/nuSQuIDS/cuda/detail/propagator.cuh` — new uploadSharedData signature
- `include/nuSQuIDS/cuda/detail/kernels.cuh` — InteractionDataGPU* in launchEvolve
- `src/cuda/cuda_backend.cu` — interaction data passthrough
- `src/cuda/propagator.cu` — full uploadInteractionData implementation
- `src/cuda/kernels.cu` — anticommutator, invlen, GammaRho, NC cascade, tau regen, Glashow, interacting RK4
- `src/nuSQuIDS.cpp` — GPU interaction data extraction
- `include/nuSQuIDS/nuSQuIDS.h` — nuSQUIDSAtm GPU interaction data extraction

**Known limitations:**
- Lagged cascade approximation (nc_factors computed at step start, not per sub-step)
- Tau regeneration uses O(ne³) triple loop (could be optimized with intermediate accumulator)
- Not yet validated against CPU — needs test/cuda/cuda_interactions_test.cpp
- SU(4)+ interaction kernels not implemented (SU(3) only)

**Next steps:**
- Create GPU interactions test (CPU vs GPU comparison with interactions=true)
- Deploy to cluster and validate
- Profile shared memory usage and register pressure
- Optimize tau regen to O(ne²) with intermediate tau flux array in shared memory
