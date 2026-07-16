# `cuda-backend` → `master` merge checklist

**Branch state:** 86 commits ahead of `origin/master`, 0 behind.
**Last tip:** `3cf8b80` (Fail-fast on SU(4-6) instead of silent no-op).
**Recommended target commit count after rebase:** 22-25.

## Pre-merge code hygiene — DONE

- [x] **`MAX_INT_PAIRS` host-side precondition + device `__trap`** (`bb138c1`).
  Latent silent-wrong-physics if `nrhos * ⌈ne/blockDim.x⌉ > 4`. Now throws at host before launch.
- [x] **Remove debug `printf` residue** (`0a845f2`). Per-launch "OK: completed in N steps" was bloating SLURM logs.
- [x] **SU(4-6) fail-fast** (`3cf8b80`). Was silent no-op when caller set `numneu>3`; now throws clear error.

## Pre-merge code hygiene — DEFERRED (not blocking)

- [ ] Multi-path-per-block scheduling (would unblock IceCube-shape workloads at ncz<108). Tracked as future work in the paper.
- [ ] Embedded RK pair on GPU integrator. Attempted as DOPRI5 in `d693ca2`; reverted in `bc28d80` after 9× regression on IceCube full-int. Future work.

## Commit history classification

86 commits total, fall into 11 logical phases. Recommended squash structure:

| Phase | Commits | Action | Resulting commit message |
|---|---|---|---|
| 1. Initial GPU backend (oscillation only) | `dbc9029`, `14feb56`, `70584b8`, `df7fe3e` | **squash 4 → 1** | Add CUDA GPU backend with oscillation-only nuSQUIDSAtm and single-baseline support |
| 2. Interactions infrastructure | `3ff6cc7`, `efe9fd8`, `e663002`, `1c097d9`, `bf5495e`, `62502e9` | **squash 6 → 1** | Add GPU non-coherent interactions (NC, CC, tau regen, Glashow) |
| 3. Launch + compilation fixes | `2a33eba`, `bb68e5c`, `512f210`, `c0d05a0`, `e053945`, `7281526`, `3f912c8`, `aabfa70` | **squash 8 → 1** | Fix GPU launch infrastructure (stack size, MIG handling, error codes) |
| 4. Lambda/HI refactors | `2200800`, `90e8d26` | **keep both** | (already named well) |
| 5. Debug-print scaffolding | `be71595`, `6937e4c`, `a3c4742`, `2a6bc75`, `6b9804a`, `a488379`, `c4c5349`, `f39cd85`, `fb33e33`, `574204a`, `f7610dc`, `b285a9a`, `de00be6`, `4db5a1c`, `4aa83b4`, `5436216`, `2fbc197` | **DROP all 17** | (transient debug, no surviving signal) |
| 6. Unit conversion (cm/eV → natural) | `b5379ab`, `39c494c`, `829ba3c`, `97da3c7`, `1bdd69a`, `3acee6b` | **squash 6 → 1** | Use SQuIDS Const natural units throughout GPU interaction kernels |
| 7. Diagnostic test iteration | `8c34d05`, `39401f7`, `494438c`, `b551211`, `ba9cbf7`, `8237bb2`, `cf111a6`, `20f73b0`, `3669f5a`, `2554360`, `c04d005` | **squash 11 → 1** | Add cuda_rhs_diagnostic test for per-component CPU/GPU comparison |
| 8. Correctness fixes (load-bearing) | `f3ec3a8` (rounded_ne), `2b19c42` (suTrace3 normalization), `f041533` (cascade refresh between half-steps), `cf6fc71` (test OOB), `1bf7012` (constant-density target ndens), `43b3d12` (Glashow units), `f3b9a29` (tau regen target ndens), `4e10c46` (single_baseline Earth track), `d8bc29f` (tau regen rho index), `72d2c14` (Glashow regression test), `efac5ea`/`1bcd2b8` (test infra) | **keep all individually** | (each is its own physics finding worth preserving) |
| 9. MAX_PAIRS evolution | `9dae151`, `27485b9`, `99fab1f` | **squash 3 → 1** | Perf #1: Move RK4 staging buffers from per-thread stack to per-block workspace |
| 10. Final correctness + Perf #2 | `1f619ec` (Richardson consistency), `655080d` (Glashow body-composition), `73eb097` (substage cascade refresh), `4596c92` (OR-pass test criterion), `f7f61b4` (E^-2/E^-3 tests), `3bab4cc` (Perf #2 profile cache) | **keep all individually** | |
| 11. Process / admin / churn | `287afa6` (Codex review doc), `2f30ea4` (subagent defs), `3a3445a` (plan), `e087758` (Codex re-review), `a1dfdef` (merge master), `fb83a51`/`7da25d3` (Perf #3 +revert), `09b2038` (paper-diag infra), `d693ca2`/`bc28d80` (DOPRI5 +revert), `bb138c1`/`0a845f2`/`3cf8b80` (this cleanup) | **mostly keep**; consider dropping `287afa6` + `e087758` (Codex review docs) and the two revert pairs (`fb83a51`+`7da25d3`, `d693ca2`+`bc28d80` — zero-sum) | |

**Target after rebase: ~22-25 commits.** Each tells a coherent physics or perf story.

## Recommended `git rebase -i` plan

```bash
# Start a rebase from master with all commits as `pick`s
git checkout cuda-backend
git rebase -i origin/master
```

Then in the editor, group the lines per the table above. Examples:

```
pick   dbc9029 Add CUDA GPU backend for nuSQUIDSAtm
squash 14feb56 Move CUDA tests to test/cuda/ and add detailed diagnostics
squash 70584b8 Match cuda_test parameters to atmospheric_osc regression test
squash df7fe3e Add single-baseline GPU support and update bindings
# (combined commit message: "Add CUDA GPU backend with oscillation-only ...")

pick   3ff6cc7 Add cross-section data upload infrastructure for GPU interactions
squash efe9fd8 Add invlen computation, anticommutator, GammaRho, and interacting RK4
squash e663002 Add NC cascade kernel with cooperative shared-memory computation
squash 1c097d9 Add tau regen, Glashow cascade, and enable GPU interaction path
squash bf5495e Update WORKLOG with full GPU interaction implementation details
squash 62502e9 Add GPU interactions validation test
# (combined: "Add GPU non-coherent interactions (NC, CC, tau regen, Glashow)")

# Drop all DEBUG-prefixed commits and other transient artifacts
drop   be71595 Add kernel debug prints
drop   6937e4c Add density/state debug prints to kernel
drop   a3c4742 DEBUG: oscillation-only in interacting RK4 to isolate issue
drop   2a6bc75 DEBUG: use rk4StepSU3 everywhere to test shared-memory preamble
drop   6b9804a DEBUG: skip interaction preamble entirely
drop   a488379 DEBUG: interactions=false in test to isolate energy range issue
drop   c4c5349 DEBUG: print derivative components in interacting step
drop   f39cd85 DEBUG: track step acceptance
drop   fb33e33 DEBUG: print target number densities
drop   574204a Fix debug print location
drop   f7610dc DEBUG: print per-flavor invlen, nc_factors, Gamma, F_int at ie=0
drop   b285a9a DEBUG: print per-rho error at ie=0
drop   de00be6 DEBUG: print sf/sh/y/st values for first pair
drop   4db5a1c DEBUG: disable absorption to isolate cascade
drop   4aa83b4 DEBUG: disable cascade too — pure oscillation in interacting path
drop   5436216 DEBUG: print nc_factors and F_int magnitudes
drop   2fbc197 DEBUG: print all (ie,rho) pairs with error>100 to find culprit

# Drop the two failed perf experiments + their reverts (zero-sum)
drop   fb83a51 Perf #3: Two-pass tau regeneration (O(ne^3) -> O(ne^2))
drop   7da25d3 Revert "Perf #3: Two-pass tau regeneration (O(ne^3) -> O(ne^2))"
drop   d693ca2 Replace RK4+Richardson with DOPRI5 5(4) FSAL on the interaction path
drop   bc28d80 Revert DOPRI5 — production regression on stiff workloads

# Drop transient Codex review documents (the value is in the resulting fixes)
drop   287afa6 Add Codex code review of GPU interaction backend
drop   e087758 Add Codex re-review of cuda-backend after Perf #1

# Keep the rest as pick
```

## Pre-merge verification (run after rebase)

```bash
# Smoke test build
./configure --enable-cuda && make -j4

# Full cluster verification
bash deploy_cuda.sh
ssh mimo "tail -200 /n/holylfs05/LABS/.../test/cuda/cuda_test_*.out | grep -E 'PASS|FAIL'"

# Expected: all 5 cuda_interactions_test cases PASS, all 5 single_baseline cases PASS,
# benchmark sweep (5,20)→(200,200) completes without crashes,
# cuda_paper_diagnostics IceCube workload runs (~22-25s GPU on full A100).
```

## Post-merge

- [ ] Tag the merge commit (e.g., `v1.14.0-gpu`).
- [ ] Update `CHANGELOG.md` with the GPU backend section.
- [ ] Bump version in `configure.ac` if the project uses one.
- [ ] Run `cd test && bash run_tests` for CPU regression (the upstream test workflow added in PR #54 will run automatically on GitHub).
