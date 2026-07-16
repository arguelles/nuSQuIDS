# Codex Re-Review of cuda-backend after Perf #1

Review date: 2026-04-27
Branch: `cuda-backend` (HEAD `99fab1f`)
Scope: GPU interaction backend (correctness vs CPU oracle, performance, numerics, branch hygiene).
Prior review (stale): `CODEX_REVIEW.md`, commit `287afa6`. Findings already closed by `f041533`, `1bf7012`, `43b3d12`, `2b19c42` are NOT re-flagged here unless they regressed.

Note on methodology: the OpenAI Codex MCP endpoint was unavailable for this session (permission denied). The review below was conducted by manual file-by-file inspection of every file in scope plus the CPU oracle in `src/nuSQuIDS.cpp`. Findings cite verified file:line locations.

---

## Summary

| Severity | Count |
|---|---|
| blocker | 1 |
| important | 6 |
| nit | 11 (5 of which are explicit "no issues found" entries documenting category coverage) |

All four review categories produced findings. The single blocker is a **silent SU(4) no-op** that would let a user request `numneu=4` and get a kernel that immediately returns without any error or evolution. Two important correctness items were also identified: (a) Glashow electron-density formula diverges from CPU when `body_has_composition`, and (b) cascade refresh between the full step and the half steps still uses different `nc_factors` (the Richardson estimator compares two different ODEs). The remaining important items are performance/hygiene risks for the upcoming master merge.

---

## Correctness vs CPU oracle

### [blocker] Silent SU(4) no-op — src/cuda/kernels.cu:861, 1245
`evolveKernelImpl<NFLV>` returns immediately when `NFLV != 3` (line 861), but `launchEvolve` still dispatches `evolveKernelImpl<4>` for `numneu == 4` (line 1246). A user requesting `numneu=4` gets a kernel that copies initial states unchanged, with zero diagnostic output. Oscillation-only and interacting paths are both affected.
**Recommendation:** In `launchEvolve` (or higher up in `cuda_backend.cu::Evolve`), reject `numneu != 3` with an explicit `throw std::runtime_error("CUDABackend currently supports only numneu=3")`, and remove the `case 4:` dispatch until the SU(4) kernels are implemented.
**Route to:** physics-dev

### [important] Glashow `num_e` ignores body composition path — src/cuda/kernels.cu:239-246, 633-639
`electronNumberDensitySU3` only handles the isoscalar (`n_targets==1 → ndens[0]*ye`) and pure p/n (`n_targets>=2 → ndens[0]`) cases. CPU `nuSQUIDS::UpdateInteractions` lines 772-784 takes a different branch when `body_has_composition || hasNuclearXS` and computes `num_e = density / (m_e + m_p + m_n*(1-ye)/ye)`, derived from mass density. For nuclear-target cross sections, `target_ndens[0]` on GPU is the first nuclear target's number density (e.g. ¹⁶O), not the proton/electron density. The two formulas only agree for plain p/n media without composition.
**Recommendation:** Either (a) plumb the CPU-computed `num_e` through `GPUPathData` as a per-sample spline alongside `target_ndens`, then evaluate it directly on GPU, or (b) replicate the `body_has_composition` and `hasNuclearXS` branches on GPU using auxiliary inputs. Prefer (a) — the user has explicitly mandated using `GetTargetNumberDensities()` and friends on the CPU side and shipping the result to GPU.
**Route to:** physics-dev

### [important] Richardson estimator compares two different ODEs — src/cuda/kernels.cu:1027-1124
The interacting RK4 path computes `sf` (full step h) using `nc_factors` snapshotted at `x` (lines 1043-1046), then refreshes `nc_factors` at `x_mid` and computes `sh` (two half-steps, second half at `x+h/2`) using the refreshed `nc_factors` (line 1108-1111). The error estimate `|sf - sh| / 15` and the corrected state `sh + (sh - sf)/15` therefore mix two different approximate ODEs, not two different step sizes of the same ODE. Richardson extrapolation assumes a common ODE. The prior Codex review flagged this and the cascade refresh between half-steps (`f041533`) only partially addresses it — `sf` is still computed with start-of-step factors.
**Recommendation:** Either (a) recompute `sf` after the midpoint refresh so both estimates use the refreshed factors as the reference (cheap, but loses the "compare two step sizes" semantics for the lagged-cascade approximation), or (b) make cascade refresh a per-RK-stage operation so all four stages of both `sf` and `sh` see the cascade evaluated at the matching `(x, y)` pair (CPU-equivalent, expensive). Option (a) restores Richardson validity at the cost of one more derivative pass; option (b) is the proper fix and aligns with the `PreDerive`-on-every-RHS CPU semantics.
**Route to:** physics-dev

### [nit] Tau regen recomputes production sum redundantly per (e1, et) — src/cuda/kernels.cu:518-574
GPU's `computeTauRegenSU3` puts the production loop `for (int en = et+1; en < ne; en++)` *inside* the `for (int et = e1+1; et < ne; et++)` loop (line 527), recomputing `tau_decay_fluxes[et]` once per (e1, et) pair instead of once per et. CPU computes it once and stores into `tau_decay_fluxes[et]` then iterates et to add to `tau_lep_decays[e1]`. The mathematical result is identical (verified by tracing both formulas against `src/nuSQuIDS.cpp:1053-1106`), but the GPU does O(ne) redundant work per cascade refresh. Not a correctness bug; flagging as a performance concern carried over to perf-dev once Task 3's residual is closed.
**Recommendation:** After Task 3 closes the 40% tau regen residual, refactor to pass 1 (compute `tau_decay_fluxes[et]` for all et into shared memory) followed by pass 2 (accumulate decays into `nc_factors[e1]`). Not load-bearing for correctness — just a perf nit until then.
**Route to:** perf-dev

### [nit] No issues in `rounded_ne` vs `ne` strides — src/cuda/kernels.cu, interactions_gpu.cuh
Verified each cross-section access against the CPU layout:
- `sigma_NC/CC` shape `[trg][rho][flv][rounded_ne]` → GPU `sigma_index` last dim = `rounded_ne` (interactions_gpu.cuh:84).
- `dNdE_NC/CC` shape `[trg][rho][flv][ne][rounded_ne]` → GPU `dNdE_index` middle dim = `ne`, last dim = `rounded_ne` (interactions_gpu.cuh:77).
- `dNdE_tau_all/lep` shape `[nrhos][ne][rounded_ne]` → GPU uses `(rho*ne+et)*rne+e1` (kernels.cu:586-587).
- `dNdE_GR` shape `[ne][rounded_ne]` → GPU uses `e2*rounded_ne+e1` (kernels.cu:676).
- `sigma_GR` shape `[rounded_ne]` (1D) → GPU uses `[e2]` directly (kernels.cu:673), upload size correct (`propagator.cu:137`).
All consistent. No change needed.
**Recommendation:** None.
**Route to:** none

### [nit] No `__device__` lambda capture risk in cascade/derivative path
Verified: no `[=]` or `[&]` device lambdas in `kernels.cu`. The previous review's concern about nvcc silently zeroing lambda-captured state was addressed by the explicit `computeDerivativeSU3` device function (line 717 comment confirms this).
**Recommendation:** None.
**Route to:** none

### [nit] `convertPath` constant-density path correctly populated — src/cuda/cuda_backend.cu:147-166
Verified: `profile.n_targets = path.n_targets` and `constant_target_fractions` is filled with `path.target_ndens[t][0]` (number density in eV³) on the CONSTANT branch. `evaluateTargetFraction` returns `constant_target_fractions[t]` for non-TABULATED profiles (`body_gpu.cuh:162-167`), so `computeNCCascadeSU3` and `computeInvlenSU3` see the correct ndens. The fix from `1bf7012` is intact.
**Recommendation:** None.
**Route to:** none

---

## Performance

### [important] Per-thread stack still ~85 doubles in interacting path — src/cuda/kernels.cu:776-829, 717-765
Even after Perf #1 moved the `corrected_buf`/`sf_buf` to global workspaces, the inner RK4 step still allocates per thread:
- In `evolveKernelImpl` interacting path (line 1037-1063): `y[9]`, `st[9]`, `sf[SU]` is a workspace pointer (good), `sh[9]` (line 1107). ≈27 doubles.
- In `rk4StepSU3_interacting` (line 788): `k[9], acc[9], tmp[9]` = 27 doubles.
- In `computeDerivativeSU3` (line 731-756): `HI[9]`, `target_ndens[16]`, `invlen[3]`, `evol_proj[27]`, `Gamma[9]`, `acomm[9]`, `F_int[9]` = 82 doubles.
Total ≈ 136 doubles × 8 B = ~1.1 KB per thread *just for the derivative call frame*. With 128 threads × 1.1 KB = 140 KB per block of stack frame at peak depth (excluding any nested `iCommutatorSU3` / `antiCommutatorSU3` temporaries which are a few more 9-component locals). The user-set `cudaLimitStackSize=8192` (`propagator.cu:190`) is well above this, but local-memory traffic from these big frames likely contributes most of the residual register pressure that Perf #1 was meant to address.
**Recommendation:** The biggest single contributor is `target_ndens[MAX_TARGETS=16]` allocated per-thread inside `computeDerivativeSU3` and recomputed per RK substep. Hoist target ndens into shared memory once per macro-step (or per cascade refresh), and read it from shared in the derivative path. Second largest is `evol_proj[27]` — also a candidate for shared since H0 and projectors don't depend on the substep state.
**Route to:** perf-dev

### [important] Cascade kernels do redundant target spline / projector evaluation — src/cuda/kernels.cu:398-406, 502-510, 635-637, 738-739
Each of `computeNCCascadeSU3`, `computeTauRegenSU3`, `computeGlashowCascadeSU3`, AND every call to `computeDerivativeSU3` (4 per RK4 step × 2 half-steps) independently re-evaluates the target-density Akima splines via `evaluateTargetFraction(profile, t, x_eval)`. For a 200-energy macro-step on a body with 4 targets, that's `4 targets × (3 cascades + 4*2 derivative calls) = 44` redundant spline evaluations per energy node per macro-step. Same pattern for `evol_proj` recomputation in the derivative.
**Recommendation:** Batch the spline evaluations: compute `target_ndens[t]` and `evol_proj[flv][9]` once per cascade refresh into shared/global memory, and have all consumers read from there. This is the natural fit for the planned Perf #2 (Task 4) and dovetails with the Perf #3 `src_factor` workspace (Task 5).
**Route to:** perf-dev

### [nit] `__launch_bounds__(128)` chosen without explicit minBlocksPerSM — src/cuda/kernels.cu:848
`__launch_bounds__(128)` allows the compiler to assume max 128 threads per block but does not constrain occupancy. With ~85 per-thread doubles in the deep frame, registers/thread is likely > 100 (the 255-reg cap is very plausible per the prior review's PTXAS output) and the kernel is occupancy-limited. A tuple `__launch_bounds__(128, 2)` or dropping to 64 threads/block could free registers.
**Recommendation:** Once Task 1 publishes the PTXAS register count, perf-dev should sweep `(threads, minBlocks) ∈ {(128,1), (128,2), (64,2), (64,4)}` and pick the best.
**Route to:** perf-dev

### [nit] Workspace allocations not reused across `Evolve()` calls — src/cuda/propagator.cu:264-270
The workspace grows on demand (`if state_bytes > workspace_size_bytes_`), but in the multi-batch single-GPU path (line 256-270), each call to `evolveBatch` re-runs `allocateBatch`. The `cudaMalloc` for `d_states_` and `d_paths_` is freed and reallocated each batch (lines 256-260), which is fine since their size depends on `n_paths`. The workspace is size-checked properly, but for a sequence of equally-sized batches there's still no issue. Just confirming this is intentional.
**Recommendation:** None — current behavior is correct. Could improve by caching `d_states_` size too if that becomes a hotspot.
**Route to:** none

---

## Numerics

### [nit] No literal cm/eV/MeV unit constants found in CUDA code
Searched all `.cu`/`.cuh` files for `cm`, `1e3`, `1e6`, `1e9`, `MeV`, `GeV`. The only matches are: numerical epsilons (`1e-30` for vacuum-density check at `cuda_backend.cu:128`, `1e-10` for ye threshold at `kernels.cu:181`, `1e-15` for time-step end check at `kernels.cu:938`, `1e-10` for `solver_config.h_min` at `cuda_backend.cu:231`), all of which are dimensionless tolerances. The `HI_constants = sqrt2*GF*Na/cm³` is computed on CPU and plumbed via `params.HI_constants` (`nuSQuIDS.cpp:134, 252, 2555`). All target number densities arrive as natural-units eV³ from `GetTargetNumberDensities()`. Unit-hygiene contract holds.
**Recommendation:** None.
**Route to:** none

### [nit] Stale `SQRT2_GF_NA = 7.63247e-14` constant — include/nuSQuIDS/cuda/detail/physics.cuh:21
Declared but never used. `grep -rn SQRT2_GF_NA` finds only the declaration. The value is also slightly stale (a 5-digit truncation of the actual SQuIDS constant; comparing against `params.sqrt2*params.GF*params.Na`). It's a foot-gun — a future contributor could pick it up thinking it's the canonical GPU value.
**Recommendation:** Delete the constant.
**Route to:** self (branch hygiene cleanup)

### [nit] No mixed fp32/fp64 boundaries
All GPU data paths verified to use `double`: state arrays, cross sections, projectors, splines, RK4 staging, tau decay spectra. No `float` literals or `__float2double` casts in the data path. The `sincos()` call in `evolveProjectorsSU3` and friends uses the double-precision overload via inferred type.
**Recommendation:** None.
**Route to:** none

### [nit] Richardson scale guard is correct — src/cuda/kernels.cu:1018-1022, 1120-1124
The error normalization `scale = abs_error + rel_error * max(|y|, |sh|)` plus the `if (scale > 0.0)` guard is well-formed. With user-set tolerances both > 0, scale is bounded below by `abs_error`. The guard prevents NaN if `abs_error == 0 && rel_error == 0`, but in that pathological case the controller silently skips error accumulation, which is the right behavior (no rejection).
**Recommendation:** None.
**Route to:** none

---

## Branch hygiene

### [important] Always-on `printf("OK: completed in %d steps...")` per kernel launch — src/cuda/kernels.cu:1171-1173
The "OK" branch fires on every successful launch (gated only by `path_idx==0 && threadIdx.x==0`, which is always true once per launch). For a benchmark that calls `Evolve()` repeatedly, this floods stderr and serializes via the GPU printf buffer, distorting wall-time measurements. The "WARNING: max_steps reached" branch is legitimate diagnostic; the "OK" branch is debug residue from validation.
**Recommendation:** Delete the `else` branch (lines 1171-1173). Keep the `WARNING` branch (a real error condition users should see).
**Route to:** self (branch hygiene cleanup)

### [important] `physics.cuh` template family declared but never defined or used — include/nuSQuIDS/cuda/detail/physics.cuh:77-168
`computeH0`, `computeHI`, `computeGammaRho`, `computeInteractionsRho`, `evolveProjectors`, `updateInteractionLengths` are template `__device__` declarations with no definitions anywhere in the codebase. The active kernels in `kernels.cu` use their own `computeHI_SU3`, `computeGammaRhoSU3`, `evolveProjectorsSU3`, `computeInvlenSU3`, `computeInteractionsRhoSU3` instead. Anyone reading `physics.cuh` will assume these templates are the architecture; finding them dead is a maintenance trap.
**Recommendation:** Either delete the unused template declarations (preferred — they were clearly an early API draft), or move them to a `_planned_` block with a TODO referencing the SU(4-6) effort.
**Route to:** self (branch hygiene cleanup)

### [important] `sumath.cuh` SU(N) commutator/anticommutator/evolve declarations — include/nuSQuIDS/cuda/detail/sumath.cuh:121-167
`suCommutator3/4`, `suAnticommutator3/4`, `suEvolve3/4` are declared but never defined. The active kernels use inline functions `iCommutatorSU3`, `antiCommutatorSU3`, `evolveProjectorsSU3` defined directly in `kernels.cu`. Same dead-declaration problem as `physics.cuh`. The `suTrace3`/`suTrace4` definitions in this header ARE used and load-bearing.
**Recommendation:** Delete the six unused declarations. Keep `suTrace3`, `suTrace4`, `suTrace<NFLV>`, `suScale`, `suAdd`, `suAddScaled`, `suCopy`, `suZero`.
**Route to:** self (branch hygiene cleanup)

### [nit] `cudaLimitStackSize = 8192` may be reducible after Perf #1 — src/cuda/propagator.cu:190
The 8 KB per-thread stack limit was set when `corrected_buf` and `sf_buf` lived on the stack. Per-thread stack is now ~1.1 KB (see Perf finding above). Reducing the limit to 4096 or 2048 would lift MIG-partition support (the comment at line 188-189 acknowledges 8 KB sometimes fails on MIG).
**Recommendation:** After Task 1's PTXAS report, perf-dev should set `cudaLimitStackSize` to `2 * (measured stack frame)` and verify on MIG.
**Route to:** perf-dev

### [nit] No printf/TODO/FIXME residue beyond the two flagged matches
`grep -nE 'printf|TODO|FIXME|XXX' src/cuda/*.cu include/nuSQuIDS/cuda/detail/*.cuh` returned only the two `printf` lines in `kernels.cu:1169` and `kernels.cu:1172` (both addressed in the "always-on printf" finding above). No `TODO`/`FIXME`/`XXX` markers anywhere.
**Recommendation:** None.
**Route to:** none

---

## Notes for downstream work

- The **40% tau regen residual** flagged for Task 3 was NOT identified by static review — the GPU and CPU formulas trace identical line-by-line. The bug is most likely a runtime ordering issue (e.g. cascade refresh schedule between half-steps), a stride/indexing edge case at `e1=0` or `et=ne-1`, or a target-fraction issue specific to the test's body composition. Recommend the diagnostic from the plan's Task 3 step 1 (per-bin component test) before further code review.
- **Glashow `num_e` fix** (this review's important #1) is unlikely to affect the existing Glashow regression test (Earth-atmosphere, no composition data), but will trigger the moment the test suite adds a body with `composition()`. Worth fixing before merge regardless.
- **Richardson estimator validity** (this review's important #2) is a long-standing semantics issue. A pragmatic interim fix is to recompute `sf` after the cascade midpoint refresh — adds one more derivative pass per accepted step but restores estimator validity.
