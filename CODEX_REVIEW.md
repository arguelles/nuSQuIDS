# Codex Code Review — nuSQuIDS GPU Interaction Backend

Review date: 2026-04-22
Scope: CUDA backend for NC/CC interactions (post-validation at 0.2%)

## Executive Summary
- The SU(3) algebra used by the active GPU path now matches SQuIDS for `SUTrace`, `iCommutator`, `ACommutator`, and SinCos projector evolution; I did not find another definite SU(3) normalization/sign bug. The header still advertises SU(4) operations without definitions, and there are no `suTrace5/6` GPU paths to review.
- Constant-density GPU paths with interactions appear broken: target number-density splines are only populated for tabulated profiles, while constant profiles return before target densities are installed, so CC/NC inverse lengths can become zero.
- Glashow and tau-regeneration paths are not CPU-equivalent. Glashow uses `density * ye` instead of the CPU electron number density and is also omitted from absorption; tau regeneration uses `ye` target fractions and raw mass density instead of the uploaded target number densities.
- The NC cascade quadrature itself matches the CPU's upper-triangular energy sum, but the RK4 integration still uses frozen cascade factors inside each RK stage. The midpoint refresh improves the two-half-step path, but the full step and per-stage RHS calls remain inconsistent with CPU `PreDerive()`.
- The reported register/stack pressure is expected from per-thread O(ne) result buffers plus repeated 9-component temporaries. The largest stack consumers are the local RK buffers, not just the cascade math.

## Findings

### [CORRECTNESS-BUG] Constant-density interaction paths zero target densities — src/cuda/cuda_backend.cu:147
**Problem:** Constant-density GPU profiles drop the precomputed target number densities. `convertPath()` returns a `CONSTANT` profile with `profile.n_targets = 0` before fitting or storing `path.target_ndens`, while interaction kernels still loop over `interaction_data.n_targets`.

**Evidence:** The constant branch returns at `src/cuda/cuda_backend.cu:147-152`. Target-density splines are only populated later in the tabulated branch at `src/cuda/cuda_backend.cu:163-168`. During upload, `constant_target_fractions` defaults to empty/zero-backed values at `src/cuda/propagator.cu:305-308`, and `evaluateTargetFraction()` returns `profile.constant_target_fractions[target]` for non-tabulated profiles at `include/nuSQuIDS/cuda/detail/body_gpu.cuh:162-165`. The interaction code then treats these values as actual target number densities in `computeDerivativeSU3()` and `computeNCCascadeSU3()` at `src/cuda/kernels.cu:676-681` and `src/cuda/kernels.cu:365-373`.

**Proposed fix:** Do not classify a path as `CONSTANT` without also preserving target number densities when interactions are enabled. Either populate `profile.constant_target_fractions[t] = path.target_ndens[t][0]` and `profile.n_targets = path.n_targets` before the constant return, or force interacting paths through the tabulated profile path so `target_a*` splines exist.

**Risk if ignored:** CC/NC absorption and NC cascade can be zero on exactly the common constant-density test case, giving falsely excellent oscillation-only behavior and badly wrong interacting propagation.

### [CORRECTNESS-BUG] Glashow absorption/source uses the wrong density and misses `invlen_INT` — src/cuda/kernels.cu:613
**Problem:** GPU Glashow source uses `sigma_GR * density * ye` with `density` in g/cm^3, while CPU computes electron number density in natural units and stores it in `int_state.invlen_GR`. GPU absorption also never adds Glashow to electron-antineutrino `invlen_INT`.

**Evidence:** CPU builds `num_e` using target densities or mass-density conversion at `src/nuSQuIDS.cpp:754-803`, then adds `invlen_GR` into `invlen_INT` for the electron antineutrino at `src/nuSQuIDS.cpp:794-803`. GPU source uses `idata.d_sigma_GR[e2] * density * ye` at `src/cuda/kernels.cu:612-616`. GPU absorption `computeInvlenSU3()` only sums `sigma_CC + sigma_NC` at `src/cuda/kernels.cu:245-254`, and `computeGammaRhoSU3()` consumes only that result at `src/cuda/kernels.cu:680-687`.

**Proposed fix:** Add a device helper that reproduces the CPU `num_e` logic from the uploaded target number-density data/profile, then compute `invlen_GR = sigma_GR * num_e`. Add it into the electron flavor `invlen` for `rho == antineutrino` when Glashow is enabled, and reuse the same `invlen_GR` in `computeGlashowCascadeSU3()`.

**Risk if ignored:** Electron-antineutrino attenuation and Glashow regeneration are dimensionally wrong and can be orders of magnitude off outside NC+CC-only validation.

### [CORRECTNESS-BUG] Tau regeneration uses mass density and `ye` fractions instead of CPU target number densities — src/cuda/kernels.cu:464
**Problem:** `computeTauRegenSU3()` approximates target fractions as `{ye, 1-ye}` and multiplies the summed cross section by raw `density`. CPU uses `int_state.invlen_CC`, which has already folded in per-target number densities from `GetTargetNumberDensities()`.

**Evidence:** GPU constructs `tfrac` from `ye` at `src/cuda/kernels.cu:464-466`, sums `tfrac[t] * sigma_CC` at `src/cuda/kernels.cu:511-515`, and then multiplies by `density` at `src/cuda/kernels.cu:516`. CPU computes `invlen_CC` from `sigma_CC * targetDensities[trg]` at `src/nuSQuIDS.cpp:740-749`, then tau regeneration uses `int_state.invlen_CC[...]` plus target fractions in the cascade weights at `src/nuSQuIDS.cpp:1058-1064`.

**Proposed fix:** Pass/evaluate the same target number densities used by `computeInvlenSU3()` into tau regeneration. For each source rho/flavor/energy, compute `invlen_CC_tau = sum_t target_ndens[t] * sigma_CC[t,rho,tau,en]`, and use separately computed CPU-equivalent target fractions only for the differential spectrum weighting.

**Risk if ignored:** Tau regeneration is not validated by the NC+CC result and will be wrong for natural-unit scaling, nuclear targets, composition-dependent media, and even proton/neutron media where mass density is not a number density.

### [CORRECTNESS-RISK] Midpoint cascade refresh is only a half-fix for CPU RHS semantics — src/cuda/kernels.cu:1007
**Problem:** CPU calls `PreDerive()` on every RHS evaluation; that evolves projectors and refreshes interactions at the current solver stage. GPU refreshes cascade factors at the step start and once at the midpoint, but each `rk4StepSU3_interacting()` still holds one `nc_factors` array constant across all four RK stages.

**Evidence:** CPU `RHS()` calls `Derive()`, which calls `PreDerive(t)` at `SQuIDS/src/SQuIDS.cpp:548-552` and `SQuIDS/src/SQuIDS.cpp:479-496`; nuSQuIDS `PreDerive()` calls `EvolveProjectors()` and `UpdateInteractions()` at `src/nuSQuIDS.cpp:278-296`. GPU `rk4StepSU3_interacting()` passes the same `nc_factors` into all four derivative evaluations at `src/cuda/kernels.cu:726-759`. The new midpoint refresh happens only between the two half steps at `src/cuda/kernels.cu:1007-1029`; the full-step `sf_buf` was already computed with start-of-step cascade factors at `src/cuda/kernels.cu:981-984`.

**Proposed fix:** For CPU-equivalent semantics, make cascade refresh a per-stage operation: compute factors for k1 at `(x, y)`, k2/k3 at `(x+h/2, y_tmp)`, and k4 at `(x+h, y_tmp)`, with block-wide synchronization between stage state writes and cascade recomputation. A cheaper intermediate option is to drop Richardson for interacting mode and use an embedded RK pair whose accepted solution and error estimate share the same staged cascade refresh schedule.

**Risk if ignored:** The accepted state can be accurate for smooth cascades but the error estimator compares two different approximate ODEs: `sf` uses start-frozen cascade, `sh` uses start/mid cascade. This can under- or over-estimate error near sharp density or flux changes.

### [CORRECTNESS-RISK] Fixed per-thread RK buffers can overflow for large `ne * nrhos / blockDim` — src/cuda/kernels.cu:857
**Problem:** Each thread stores results for at most `MAX_PAIRS = 32` energy/rho pairs, but the loops cover `ceil(ne / 128) * nrhos` pairs per thread. No guard prevents `pair_idx` from exceeding 32.

**Evidence:** `corrected_buf` is declared as `double corrected_buf[MAX_PAIRS * SU]` at `src/cuda/kernels.cu:854-858`, and interaction mode also declares `sf_buf[MAX_PAIRS * SU]` at `src/cuda/kernels.cu:961-965`. `pair_idx++` is used in loops over all `ie = threadIdx.x; ie < ne; ie += blockDim.x` and all `rho` at `src/cuda/kernels.cu:926-959`, `src/cuda/kernels.cu:968-1003`, and `src/cuda/kernels.cu:1033-1068`. With 128 threads, `nrhos=2`, and `ne > 2048`, `pair_idx` can exceed 32.

**Proposed fix:** Either enforce `ceil(ne / blockDim.x) * nrhos <= MAX_PAIRS` before launch, allocate the correction buffers in shared/global memory indexed by `(rho, ie)`, or process the per-thread energy list in chunks of `MAX_PAIRS`.

**Risk if ignored:** Silent stack overwrite in high-resolution energy grids, producing non-local numerical corruption.

### [CORRECTNESS-RISK] SU(N) GPU algebra surface is narrower than the headers imply — include/nuSQuIDS/cuda/detail/sumath.cuh:87
**Problem:** The active CUDA kernel returns immediately for `NFLV != 3`, yet headers declare SU(4) commutator, anticommutator, and evolution routines and the review request mentions SU(5/6) traces. There are no GPU `suTrace5` or `suTrace6` functions to compare.

**Evidence:** `evolveKernelImpl()` returns unless `NFLV == 3` at `src/cuda/kernels.cu:794-795`. `sumath.cuh` implements only `suTrace3()` and `suTrace4()` at `include/nuSQuIDS/cuda/detail/sumath.cuh:76-96`; it only declares `suCommutator3/4`, `suAnticommutator3/4`, and `suEvolve3/4` at `include/nuSQuIDS/cuda/detail/sumath.cuh:121-167`. SQuIDS generic trace uses `dim*a0*b0 + 2*sum` at `SQuIDS/include/SQuIDS/detail/ProxyImpl.h:262-312`; the implemented GPU SU(3/4) traces match that convention.

**Proposed fix:** Make the limitation explicit in the public backend path: reject `numneu != 3` before launch with a clear error, or implement and test SU(4/5/6) operations against SQuIDS generated files. Remove or define unused declarations to avoid accidental linkage surprises.

**Risk if ignored:** Users can request unsupported dimensions and get a no-op kernel rather than a correct CPU-equivalent evolution.

### [PERF-HIGH] Local RK buffers dominate stack spill and occupancy loss — src/cuda/kernels.cu:854
**Problem:** The reported 255 registers and stack spill are consistent with large per-thread local arrays. Interaction mode has `corrected_buf` plus `sf_buf`, and derivative/cascade helpers add many 9-component temporaries.

**Evidence:** `corrected_buf[MAX_PAIRS * SU]` is 288 doubles per thread at `src/cuda/kernels.cu:854-858`; interaction mode adds another 288-double `sf_buf` at `src/cuda/kernels.cu:961-965`. Each RK call uses `k[9], acc[9], tmp[9]` at `src/cuda/kernels.cu:724`, and each derivative adds `HI[9]`, `evol_proj[27]`, `Gamma[9]`, `acomm[9]`, and `F_int[9]` at `src/cuda/kernels.cu:669-701`. Cascade helpers add `target_ndens[MAX_TARGETS]`, `target_frac[MAX_TARGETS]`, and `evol_proj_alpha[9]` at `src/cuda/kernels.cu:365-393`.

**Proposed fix:** Move accepted/candidate state buffers to shared memory or a temporary global workspace indexed by `(path, rho, ie, component)`, and split the interacting solve into cooperative stage kernels: stage-state write, cascade factor kernel, derivative/update kernel. Also consider `__launch_bounds__(64, minBlocksPerSM)` after reducing local arrays, and mark cold features (tau/Glashow) as separate kernels to keep the NC+CC kernel smaller.

**Risk if ignored:** Low occupancy and local-memory traffic will erase much of the GPU advantage, especially with interactions enabled.

### [PERF-LOW] Shared memory layout is reasonable, but state stride and access order are not ideal — src/cuda/kernels.cu:860
**Problem:** The shared layout is compact and `s_nc_factors[(rho*3+flv)*ne + ie]` gives contiguous stores across `ie`, but `s_state[(rho*ne+ie)*9+c]` means each SU vector starts 9 doubles after the previous one. The cascade also recomputes evolved projectors and scans upper-triangular energies per output energy, causing divergent per-thread work.

**Evidence:** Shared memory layout is declared at `src/cuda/kernels.cu:860-867`; NC factor stores are contiguous over `e1` at `src/cuda/kernels.cu:439`, and reads are strided by `ne` over flavor at `src/cuda/kernels.cu:640-645`. State loads into shared are linearly coalesced at `src/cuda/kernels.cu:880-884`, but cascade reads `state_e2 = s_state + (rho * ne + e2) * SU` inside a triangular loop at `src/cuda/kernels.cu:386-418`.

**Proposed fix:** Keep the NC factor rank ordering for coalesced writes, but benchmark padding SU(3) state vectors to 10 or 16 doubles in shared for alignment and simpler vectorized loads. A larger improvement is to invert the NC cascade loop: assign cooperative work over `(rho, alpha, e2)` and accumulate into `e1 < e2`, using block reductions or atomics, so evolved projectors and `flux * invlen * dE` are computed once per source energy.

**Risk if ignored:** This is mostly performance risk; the current layout is not an obvious correctness problem and bank conflicts are likely secondary to register/local-memory pressure.

### [CORRECTNESS-RISK] RK4 step-doubling is not the CPU solver and may be weak near stiff features — src/cuda/kernels.cu:774
**Problem:** CPU SQuIDS defaults to GSL RKF45, while GPU uses classical RK4 step doubling with Richardson extrapolation. RK4 step doubling is valid for smooth non-stiff problems, but it is not bit-for-close to the CPU method and provides a less informative embedded error estimate for interaction terms that depend nonlinearly on the full state.

**Evidence:** CPU default `step` is `gsl_odeiv2_step_rkf45` at `SQuIDS/src/SQuIDS.cpp:45`, and adaptive evolution uses `gsl_odeiv2_driver_apply()` at `SQuIDS/src/SQuIDS.cpp:509-529`. GPU documents and implements full-step/two-half-step RK4 with Richardson at `src/cuda/kernels.cu:774-779` and `src/cuda/kernels.cu:948-956`/`src/cuda/kernels.cu:1055-1065`.

**Proposed fix:** Implement an embedded RK45 pair for the GPU interaction path. Dormand-Prince 5(4) costs 7 RHS evaluations per accepted step; Cash-Karp costs 6. The main cost is not the algebra but staging global/shared state and recomputing cascade factors per stage. Expect roughly 6-7 cooperative cascade refreshes per attempted step for CPU semantics; if that is too expensive, provide `solver=cpu_equivalent_rkf45` and `solver=fast_rk4_richardson` modes and validate the latter over resonance-heavy cases.

**Risk if ignored:** The GPU can pass broad NC+CC regression tests while missing localized resonance or sharp-density errors that RKF45 catches differently.

### [PERF-LOW] `suTrace3` precision is probably adequate, but deterministic compensated sums should be tested — include/nuSQuIDS/cuda/detail/sumath.cuh:79
**Problem:** `suTrace3()` is a short 9-term double sum. Kahan or pairwise summation could reduce cancellation when extracting tiny flavor components, but it will not fix a physics-level discrepancy and may increase register pressure.

**Evidence:** GPU `suTrace3()` uses `double result = 3.0*a[0]*b[0]` plus eight sequential terms at `include/nuSQuIDS/cuda/detail/sumath.cuh:79-84`. SQuIDS generic trace uses the same mathematical convention but a loop order that depends on odd/even dimension handling at `SQuIDS/include/SQuIDS/detail/ProxyImpl.h:262-312`. The reviewed CUDA files use `double` state, cross-section, profile, and solver values; I did not find `float` in the reviewed backend data path.

**Proposed fix:** Leave default `suTrace3()` as-is for performance, but add a compile-time option or test helper using pairwise summation:
`3*a0*b0 + ((2*a1*b1 + 2*a2*b2) + ...)`.
Use it in CPU/GPU comparison tests to quantify whether trace extraction contributes to the 0.2% residual. Do not enable Kahan by default until register impact is measured.

**Risk if ignored:** Low. Possible last-digit differences remain, but the current implementation is mathematically consistent with SQuIDS.

### [STYLE] Device headers contain stale comments and unused abstractions — include/nuSQuIDS/cuda/detail/physics.cuh:105
**Problem:** Several headers describe APIs or signs that the actual implementation does not use, making future fixes more error-prone.

**Evidence:** `computeGammaRho()` is documented as `GammaRho = -0.5 * ...` at `include/nuSQuIDS/cuda/detail/physics.cuh:105-107`, while the actual SU(3) helper correctly builds positive `Gamma` and subtracts the anticommutator later at `src/cuda/kernels.cu:267-280` and `src/cuda/kernels.cu:698-701`. `physics.cuh` declares generic `computeInteractionsRho()` and `updateInteractionLengths()` at `include/nuSQuIDS/cuda/detail/physics.cuh:120-168`, but the reviewed implementation uses bespoke SU(3) functions in `kernels.cu`.

**Proposed fix:** Align comments with implemented signs and either remove stale generic declarations or add definitions/tests. Keep the source of truth close to the kernel implementation.

**Risk if ignored:** Maintainers can reintroduce sign or normalization mistakes while following misleading header comments.

## Priority Order for Fix Agents
1. Fix target number-density handling for constant-density interaction profiles and add a constant-density NC+CC regression test.
2. Fix Glashow absorption/source units and `invlen_INT` inclusion; add electron-antineutrino Glashow tests against CPU.
3. Fix tau regeneration to use uploaded target number densities and CPU-equivalent target fractions; add tau-regeneration CPU/GPU tests.
4. Decide the intended interacting solver semantics: CPU-equivalent per-stage cascade refresh with RKF45/RK45, or documented fast approximation with validation bounds.
5. Move large RK local buffers out of per-thread stack and split cold cascade features to reduce register pressure.
6. Add guardrails for unsupported `numneu != 3` and for `MAX_PAIRS` overflow.

## Positive Observations
- The active SU(3) `iCommutatorSU3()` and `antiCommutatorSU3()` match the SQuIDS generated SU(3) formulas I checked against `iConmutatorSU3.txt` and `AnticonmutatorSU3.txt`.
- The fixed `suTrace3()` and `suTrace4()` use the SQuIDS trace convention `dim*a0*b0 + 2*sum_{i>0}`.
- Main sigma and dNdE indexing now consistently uses `rounded_ne` through `sigma_index()` and `dNdE_index()`.
- The NC cascade energy discretization matches the CPU upper-triangular sum over `e2 > e1` with `delE[e2-1]`.
- The shared `s_nc_factors[(rho*3+flv)*ne + ie]` layout is a sensible fit for contiguous per-energy writes and the subsequent per-thread source assembly.
