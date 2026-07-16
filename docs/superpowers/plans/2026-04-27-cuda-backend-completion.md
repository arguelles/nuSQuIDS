# CUDA Backend Completion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers-extended-cc:subagent-driven-development (recommended) or superpowers-extended-cc:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Each task names an **Executor** — dispatch the matching subagent (`physics-dev`, `perf-dev`, `codex-reviewer`, or `deployer`) defined under `.claude/agents/`.

**Goal:** Close the four open work items (#8 Codex review, #9 tau regen residual, #10 perf #2/#3, #11 commit cleanup) plus the deployer benchmark for Perf #1, so `cuda-backend` is mergeable into `master` with full physics correctness and validated GPU speedup.

**Architecture:** The CUDA GPU backend is already correctness-validated for oscillation-only and NC+CC interactions. Outstanding work is (a) measuring the Perf #1 stack-buffer change, (b) closing the tau regen 40% residual that didn't move with `f3b9a29`/`d8bc29f`, (c) two queued perf optimizations, and (d) commit cleanup before merge. Each task ships independently with a real cluster verification.

**Tech Stack:** C++17, CUDA 12.4.1, SQuIDS, GSL, HDF5. Cluster: Harvard FASRC, A100 GPUs via SLURM partition `arguelles_delgado_gpu_a100`. Deploy via `bash deploy_cuda.sh`. SSH alias `mimo`. Codex MCP for code-heavy analysis (saves tokens).

---

### Task 1: Benchmark Perf #1 (RK4 staging buffers → workspace)

**Executor:** `deployer`

**Goal:** Produce a complete before/after report for commit `99fab1f` vs baseline `1bcd2b8`: register count per thread, stack bytes per thread, end-to-end CPU-vs-GPU speedup on `cuda_interactions_test` NC+CC Test 1, and confirmation that Test 1 still passes within 0.2%.

**Files:**
- Read: `src/cuda/kernels.cu`, `src/cuda/propagator.cu`, `include/nuSQuIDS/cuda/detail/propagator.cuh` (only to confirm what changed)
- Modify: none (this is a measurement task)
- Capture: `~/scratch/nusquids/perf1-before.log`, `~/scratch/nusquids/perf1-after.log`, `~/scratch/nusquids/perf1-build-{before,after}.log`

**Acceptance Criteria:**
- [ ] Both commits build with `-Xptxas -v` output captured for `evolveKernelImpl`
- [ ] Register count and stack bytes recorded for both
- [ ] Wall-clock time for `cuda_interactions_test` NC+CC Test 1 recorded for both, on full A100 (not MIG)
- [ ] Test 1 PASS confirmed on `99fab1f` with relative residual <0.2%
- [ ] Report posted using the structured format in `.claude/agents/deployer.md`

**Verify:** `ssh mimo "grep -E 'PASS|FAIL|residual' /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda/cuda_test_*.out | tail -20"` → shows Test 1 PASS with residual <0.002

**Steps:**

- [ ] **Step 1: Sync and check out the baseline commit (`1bcd2b8`)**

```bash
ssh mimo "cd /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS && git fetch && git checkout 1bcd2b8"
```

- [ ] **Step 2: Build baseline with PTXAS verbose, capture build log**

```bash
ssh mimo "bash -l -c '
  cd /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS
  module load gcc/12.2.0 cuda/12.4.1 openmpi/5.0.5 hdf5/1.14.6 gsl/2.8
  CXX=g++ ./configure --with-squids=/n/home06/carguelles/programs/local --with-hdf5=\$HDF5_HOME --enable-cuda
  make clean && make -j4 NVCCFLAGS_EXTRA=\"-Xptxas=-v -Xptxas=--warn-on-spills\" 2>&1 | tee ~/scratch/nusquids/perf1-build-before.log
'"
```

- [ ] **Step 3: Extract baseline register/stack from build log**

```bash
ssh mimo "grep -E 'evolveKernelImpl|registers|stack frame|spill' ~/scratch/nusquids/perf1-build-before.log"
```
Expected: lines like `ptxas info : Used N registers, M stack frame, ...`. Record `N` and `M`.

- [ ] **Step 4: Run baseline test on full A100 and capture wall time**

```bash
ssh mimo "cd /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda && sbatch run_cuda_test.sh"
ssh mimo "squeue -u \$USER"   # poll until job completes (do not sleep-loop locally)
ssh mimo "ls -t /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda/cuda_test_*.out | head -1 | xargs cat" > /tmp/perf1-before.log
```

- [ ] **Step 5: Repeat steps 1–4 for `99fab1f`**

```bash
ssh mimo "cd /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS && git checkout 99fab1f"
# repeat build+test, capture to perf1-build-after.log and perf1-after.log
```

- [ ] **Step 6: Confirm Test 1 still passes within 0.2% on `99fab1f`**

```bash
grep -E 'NC.*Test 1|residual|PASS|FAIL' /tmp/perf1-after.log
```
Expected: `Test 1 PASS, max relative residual = 0.001x` (must be `<0.002`).

- [ ] **Step 7: Post structured report**

Use the format from `.claude/agents/deployer.md`:

```
Commit:    99fab1f  (vs baseline 1bcd2b8)
Partition: full A100
Build:     OK
Registers: <after> / thread   (was <before>)
Stack:     <after> B / thread (was <before> B)
Test:      cuda_interactions_test NC+CC Test 1
Result:    PASS within 0.2%
Speedup:   GPU <after>× CPU   (was <before>×)
Logs:      ~/scratch/nusquids/perf1-{before,after}.log
```

- [ ] **Step 8: Hand findings back to perf-dev for Perf #2 baseline anchoring**

No commit needed; this task produces measurements, not code.

---

### Task 2: Codex review of GPU interaction backend (#8)

**Executor:** `codex-reviewer`

**Goal:** Run a fresh-eyes Codex audit of the interaction backend (NC, tau regen, Glashow) and produce a severity-ranked findings list. The previous review (`CODEX_REVIEW.md`, commit `287afa6`) is now stale — several of its findings have been fixed (`f041533`, `1bf7012`, `43b3d12`, `2b19c42`) and Perf #1 (`99fab1f`) introduced new code paths.

**Files:**
- Read: `src/cuda/kernels.cu`, `src/cuda/propagator.cu`, `src/cuda/cuda_backend.cu`, `include/nuSQuIDS/cuda/detail/sumath.cuh`, `include/nuSQuIDS/cuda/detail/interactions_gpu.cuh`, `include/nuSQuIDS/cuda/detail/kernels.cuh`, `include/nuSQuIDS/cuda/detail/propagator.cuh`
- Reference: `src/nuSQuIDS.cpp` (CPU oracle), `CODEX_REVIEW.md` (prior review)
- Create: `docs/superpowers/reviews/2026-04-27-cuda-backend-codex-review.md`

**Acceptance Criteria:**
- [ ] Each of the four review categories (correctness / performance / numerics / branch hygiene) has at least one finding or an explicit "no issues found" note
- [ ] Each finding lists: severity (`blocker` | `important` | `nit`), file:line, one-sentence description, one-sentence recommendation
- [ ] All `blocker`-severity findings are routed via plan tasks to `physics-dev` or `perf-dev`
- [ ] No raw Codex transcript pasted; only distilled findings

**Verify:** `test -f docs/superpowers/reviews/2026-04-27-cuda-backend-codex-review.md && grep -cE '^### \[(blocker|important|nit)\]' docs/superpowers/reviews/2026-04-27-cuda-backend-codex-review.md` → integer ≥ 4

**Steps:**

- [ ] **Step 1: Read all GPU backend files locally for context**

```bash
wc -l src/cuda/*.cu include/nuSQuIDS/cuda/detail/*.cuh
```

- [ ] **Step 2: Call Codex with the correctness-audit question**

Use `mcp__codex__codex` with prompt:
> Review the following nuSQuIDS GPU files for correctness against the CPU oracle in `src/nuSQuIDS.cpp`. Specifically check: (a) `rounded_ne` vs `ne` strides in cross-section indexing, (b) `suTrace3/4` normalization (should be `dim*a[0]*b[0] + 2*sum`), (c) per-step cascade refresh between RK4 half-steps, (d) lambda-capture bugs in nvcc, (e) constant-density `convertPath()` populates `constant_target_fractions`.

- [ ] **Step 3: Call Codex with the performance-audit question**

> Review for: register pressure (>200/thread is a smell), stack spills, occupancy, kernel fusion candidates between `computeNCCascadeSU3` / `computeTauRegenSU3` / `computeGlashowCascadeSU3`, redundant memory traffic across RK4 substeps. Note that Perf #1 (commit `99fab1f`) just moved RK4 staging buffers to a Propagator-owned workspace — flag any remaining per-thread stack bloat.

- [ ] **Step 4: Call Codex with the numerics-audit question**

> Review for: catastrophic cancellation, mixed fp32/fp64 paths, hardcoded conversion constants (the project uses SQuIDS `Const` natural units everywhere — flag any cm/eV/etc. literal constants in CUDA code).

- [ ] **Step 5: Call Codex with the branch-hygiene question**

> List debug `printf`s in kernels, commented-out experiments, dead SU(4-6) stub paths that could be deleted, and squash candidates among recent commits (run on commit history `git log cuda-backend ^master`).

- [ ] **Step 6: Distill into a single review document**

Create `docs/superpowers/reviews/2026-04-27-cuda-backend-codex-review.md` with one heading per finding:

```markdown
### [blocker] suTrace3 normalization off in foo.cuh:42
**Recommendation:** ...

### [important] register count 220/thread in evolveKernelImpl
**Recommendation:** ...
```

- [ ] **Step 7: Commit the review**

```bash
git add docs/superpowers/reviews/2026-04-27-cuda-backend-codex-review.md
git commit -m "Add Codex re-review of cuda-backend after Perf #1"
```

---

### Task 3: Diagnose and close tau regen 40% residual (#9)

**Executor:** `physics-dev`

**Goal:** Identify the structural bug causing the 40% residual in `cuda_interactions_test` Test 2 (tau regen). The bug survived the target-density fix (`f3b9a29`) and the dNdE rho-index fix (`d8bc29f`), so it is not a data-prep issue. Hypothesis space includes wrong rho-index mapping for tau secondary fluxes, missing `nutau` species in the source iteration, or wrong cross-section array selected at indexing time.

**Files:**
- Read: `src/cuda/kernels.cu` (`computeTauRegenSU3`), `src/nuSQuIDS.cpp` (`GetTauRegenSrc` and tau-related accessors), `include/nuSQuIDS/cuda/detail/interactions_gpu.cuh`
- Modify (likely): `src/cuda/kernels.cu`
- Modify (possibly): `include/nuSQuIDS/cuda/detail/interactions_gpu.cuh`
- Test: `test/cuda/cuda_rhs_diagnostic.cpp`, `test/cuda/cuda_interactions_test.cpp`

**Acceptance Criteria:**
- [ ] Component-level diagnostic shows GPU `computeTauRegenSU3` matches CPU `GetTauRegenSrc` to <1% per energy bin per flavor for at least 3 sampled energies
- [ ] `cuda_interactions_test` Test 2 (tau regen) PASSES with relative residual <0.2%
- [ ] Other interaction tests (NC+CC Test 1, Glashow Test 3) still PASS within 0.2% (no regression)
- [ ] If a structural change was needed (rho-index or species iteration), explain *why* in the commit message — not just *what*

**Verify:** `ssh mimo "grep -E 'Test 2.*Tau|residual' /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda/cuda_test_*.out | tail -10"` → Test 2 PASS, residual <0.002

**Steps:**

- [ ] **Step 1: Add a tighter assertion to the diagnostic harness for tau regen**

In `test/cuda/cuda_rhs_diagnostic.cpp`, isolate the tau regen contribution to the RHS for a fixed input state and compare CPU vs GPU per energy bin. Print max relative residual.

- [ ] **Step 2: Build and run the diagnostic on the cluster (via deployer)**

Hand off to `deployer` (`run_cuda_test.sh` already builds and runs `cuda_rhs_diagnostic`). Capture the per-bin residuals.

- [ ] **Step 3: If residual is per-bin uniform → suspect a global factor (units, ndens, normalization)**

Re-check `computeTauRegenSU3`:
- Target ndens propagation (already fixed in `f3b9a29`, but re-verify the access pattern matches `GetTargetNumberDensities()` shape)
- Cross-section array stride (`rounded_ne` vs `ne`, already audited but re-check for tau specifically)
- `dNdE` integration measure (`dE` vs `dlogE`)

- [ ] **Step 4: If residual is per-bin variable → suspect a structural bug (rho-index, species iteration)**

Look at how the CPU iterates rho indices in `GetTauRegenSrc`:

```bash
grep -nE "TauRegen|tau_lep|nutau" src/nuSQuIDS.cpp | head -40
```

Compare to the GPU loop in `computeTauRegenSU3`. Specifically check:
- Is the source flux read with `rho=0` (nu) and `rho=1` (nubar) handled correctly?
- Does the GPU pick up the secondary `nu_e` / `nu_mu` from tau decay with the right branching ratio?

- [ ] **Step 5: Add a single-bin printf-debug for one suspect path**

Add a `printf` guarded by `if (path_idx == 0 && energy_idx == 5)` in `computeTauRegenSU3` printing both factors of the contribution. Compare with the same logged value from `src/nuSQuIDS.cpp`.

- [ ] **Step 6: Patch the structural bug**

Apply the minimal change. Do not refactor adjacent code. If two lines fix it, change two lines.

- [ ] **Step 7: Remove the printf, rebuild via deployer, run full test suite**

```bash
bash deploy_cuda.sh
# wait for SLURM completion via squeue, then:
ssh mimo "tail -200 /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda/cuda_test_*.out | grep -E 'PASS|FAIL|residual'"
```

- [ ] **Step 8: Commit with explanation**

```bash
git add src/cuda/kernels.cu  # plus any other files touched
git commit -m "Fix tau regen <root cause in 1 phrase>

<2-3 sentences: what the bug was, why prior fixes didn't catch it, what the fix does>"
```

---

### Task 4: Implement Perf #2 — fuse `computeDerivativeSU3` stages

**Executor:** `perf-dev`

**Goal:** Eliminate the global-memory roundtrips between `computeNCCascadeSU3`, `computeTauRegenSU3`, and `computeGlashowCascadeSU3` by fusing them into a single kernel invocation per RHS evaluation, or by issuing them on overlapping CUDA streams. Expected impact (per Codex review estimate): 1.3–1.8× additional speedup on top of Perf #1.

**Depends on:** Task 1 (need post-Perf-#1 measured baseline) and Task 3 (don't optimize broken code).

**Files:**
- Modify: `src/cuda/kernels.cu` (kernel fusion), `include/nuSQuIDS/cuda/detail/kernels.cuh` (launch glue), possibly `src/cuda/propagator.cu` (stream management)
- Read: `src/cuda/kernels.cu` `computeDerivativeSU3` and the three cascade kernels

**Acceptance Criteria:**
- [ ] One of: (a) fused kernel that computes NC+tau+Glashow contributions in a single launch, OR (b) three streams with overlapping execution
- [ ] All three interaction tests still PASS within 0.2% (NC+CC, tau regen, Glashow)
- [ ] Cluster benchmark shows speedup vs Task 1's `99fab1f` baseline
- [ ] No physics-correctness code changed (suTrace, commutators, cascade quadrature, target densities are owned by `physics-dev`)

**Verify:** `ssh mimo "grep -E 'PASS|FAIL|residual|GPU.*CPU' /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda/cuda_test_*.out | tail -20"` → all three tests PASS, GPU speedup higher than Task 1 after-baseline

**Steps:**

- [ ] **Step 1: Confirm Task 1 numbers are recorded (the "before" anchor for Perf #2)**

Read the deployer report from Task 1. Note the exact register count, stack bytes, and speedup on `99fab1f`.

- [ ] **Step 2: Decide fusion vs streams via Codex consultation**

Call `mcp__codex__codex` with the three cascade kernel bodies and ask: "Given the per-thread shared state and global-memory traffic pattern, is kernel fusion or stream concurrency the better approach?"

- [ ] **Step 3: Implement the chosen approach in a single commit**

Touch only `src/cuda/kernels.cu`, `include/nuSQuIDS/cuda/detail/kernels.cuh`, and (if streams) `src/cuda/propagator.cu` + `include/nuSQuIDS/cuda/detail/propagator.cuh`.

- [ ] **Step 4: Local build to catch syntax errors**

```bash
./configure --enable-cuda && make -j4 2>&1 | tail -40
```

- [ ] **Step 5: Hand off to deployer for cluster verification**

Send `deployer` a request: build, run `cuda_interactions_test` (all 3 tests), capture `-Xptxas -v` registers/stack, measure CPU-vs-GPU wall time. Format per `.claude/agents/deployer.md`.

- [ ] **Step 6: Commit only if all three tests PASS within 0.2%**

```bash
git add src/cuda/kernels.cu include/nuSQuIDS/cuda/detail/kernels.cuh
git commit -m "Perf #2: fuse cascade kernels in computeDerivativeSU3

Single kernel launch eliminates global-memory roundtrips between
NC, tau regen, and Glashow contributions. Speedup <X>× → <Y>×
(see deployer report)."
```

- [ ] **Step 7: If any test regresses, revert and route the regression to physics-dev**

```bash
git revert HEAD
```

---

### Task 5: Implement Perf #3 — precompute NC cascade `src_factor`

**Executor:** `perf-dev`

**Goal:** The NC cascade source factor is currently recomputed inside every RK4 substep, but its dependence on the macro-step state is small enough that it can be computed once per macro-step into a Propagator-owned workspace. Expected impact: 1.1–1.3× additional speedup, plus reduced register pressure in the inner loop.

**Depends on:** Task 4 (sequence perf optimizations to measure each independently).

**Files:**
- Modify: `src/cuda/kernels.cu` (`computeNCCascadeSU3`, RK4 substep call sites), `src/cuda/propagator.cu` (allocate the new workspace), `include/nuSQuIDS/cuda/detail/propagator.cuh` (workspace member declaration)
- Read: `src/cuda/kernels.cu` to confirm where `src_factor` is currently computed

**Acceptance Criteria:**
- [ ] `src_factor` is computed once per macro-step into a `Propagator`-owned workspace, not per RK4 substep
- [ ] Workspace is sized `n_paths * nrhos * ne * SU` doubles (same shape as Perf #1 workspaces)
- [ ] All three interaction tests still PASS within 0.2%
- [ ] Cluster benchmark shows additional speedup vs Task 4's after-baseline
- [ ] No new heap allocations per launch — workspace lifetime owned by `Propagator`

**Verify:** `ssh mimo "grep -E 'PASS|FAIL|residual|GPU.*CPU' /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda/cuda_test_*.out | tail -20"` → all three tests PASS, GPU speedup higher than Task 4 after-baseline

**Steps:**

- [ ] **Step 1: Locate the current `src_factor` computation site**

```bash
grep -nE "src_factor|nc_factor|cascade.*factor" src/cuda/kernels.cu
```

- [ ] **Step 2: Verify the workspace pattern from Perf #1**

Read `include/nuSQuIDS/cuda/detail/propagator.cuh` to see how `corrected_buf_workspace` and `sf_buf_workspace` are declared and grown. Mirror the pattern for `src_factor_workspace`.

- [ ] **Step 3: Add the workspace member + allocator**

```cpp
// In propagator.cuh
double* src_factor_workspace = nullptr;
size_t src_factor_workspace_size = 0;

// In propagator.cu, in the same function that allocates the Perf #1 workspaces:
size_t needed = n_paths * nrhos * ne * SU;
if (src_factor_workspace_size < needed) {
  cudaFree(src_factor_workspace);
  cudaMalloc(&src_factor_workspace, needed * sizeof(double));
  src_factor_workspace_size = needed;
}
```

- [ ] **Step 4: Compute `src_factor` once per macro-step**

In `evolveKernelImpl` (or wherever the macro-step loop lives), call a new `precomputeSrcFactorSU3` kernel that fills `src_factor_workspace`. Then have `computeNCCascadeSU3` read from the workspace instead of recomputing.

- [ ] **Step 5: Local build**

```bash
./configure --enable-cuda && make -j4 2>&1 | tail -40
```

- [ ] **Step 6: Hand off to deployer for verification + benchmark**

Same protocol as Task 4 step 5.

- [ ] **Step 7: Commit only if tests pass and speedup is real**

```bash
git add src/cuda/kernels.cu src/cuda/propagator.cu include/nuSQuIDS/cuda/detail/propagator.cuh
git commit -m "Perf #3: precompute NC cascade src_factor into workspace

src_factor was recomputed every RK4 substep; now computed once per
macro-step into a Propagator-owned workspace. Speedup <X>× → <Y>×."
```

- [ ] **Step 8: If regression, revert**

```bash
git revert HEAD
```

---

### Task 6: Clean up debug commits before merging cuda-backend → master (#11)

**Executor:** `codex-reviewer`

**Goal:** Produce a clean, reviewable commit history on `cuda-backend` so the merge into `master` is auditable. Identify squash candidates, dead code, debug artifacts, and stub paths that should be removed or merged. Do not rewrite history beyond what's recommended by the review — the user makes the final call on the rebase.

**Depends on:** Tasks 1, 2, 3, 4, 5 (cleanup happens after the work is done).

**Files:**
- Read: full `git log cuda-backend ^master` history
- Read: `src/cuda/kernels.cu` (look for stubbed SU(4-6) paths, debug `printf`)
- Create: `docs/superpowers/reviews/2026-04-27-cuda-backend-merge-checklist.md`

**Acceptance Criteria:**
- [ ] Every commit on `cuda-backend` not in `master` is classified: `keep` / `squash-with-N` / `drop`
- [ ] Stubbed SU(4-6) kernel paths are either deleted or explicitly documented as intentional no-ops with a TODO
- [ ] No debug `printf` left in any `.cu` or `.cuh` file
- [ ] Final recommended commit count is documented (ballpark target: 10–15 commits, currently ~30)
- [ ] User has a single command they can run to perform the rebase

**Verify:** `test -f docs/superpowers/reviews/2026-04-27-cuda-backend-merge-checklist.md && grep -cE '^\| [0-9a-f]{7,40} \|' docs/superpowers/reviews/2026-04-27-cuda-backend-merge-checklist.md` → integer matching `git rev-list --count cuda-backend ^master`

**Steps:**

- [ ] **Step 1: Enumerate commits to be merged**

```bash
git log --reverse --oneline cuda-backend ^master > /tmp/cuda-backend-commits.txt
wc -l /tmp/cuda-backend-commits.txt
```

- [ ] **Step 2: Build a classification table**

Create `docs/superpowers/reviews/2026-04-27-cuda-backend-merge-checklist.md` with a table:

```markdown
| SHA | Subject | Classification | Notes |
|---|---|---|---|
| 2b19c42 | Fix suTrace3/suTrace4 normalization | keep | Load-bearing correctness fix |
| 9dae151 | Guard against MAX_PAIRS overflow | squash-with-99fab1f | Superseded |
| 27485b9 | Add compile-time MAX_PAIRS sanity check | squash-with-99fab1f | Superseded |
| ... | ... | ... | ... |
```

- [ ] **Step 3: Audit for debug `printf` and stub paths via Codex**

```bash
grep -nE 'printf|TODO|FIXME|XXX' src/cuda/*.cu include/nuSQuIDS/cuda/detail/*.cuh
```

Hand the output to `mcp__codex__codex` with: "Classify each match as (a) intentional logging that should stay, (b) debug residue to remove, or (c) a real TODO that should be tracked separately."

- [ ] **Step 4: Generate the rebase plan**

For squashes, generate the explicit `git rebase -i` instructions, e.g.:

```
pick 2b19c42 Fix suTrace3/suTrace4 normalization
pick f041533 Refresh nc_factors between half-steps
squash 9dae151 Guard against MAX_PAIRS overflow
squash 27485b9 Add compile-time MAX_PAIRS sanity check
squash 99fab1f Move RK4 staging buffers to workspace
...
```

Save this as a code block in the checklist doc.

- [ ] **Step 5: Run the post-rebase verification plan**

Document — but do NOT execute — the verification:

```bash
# After user runs the rebase:
git log --oneline master..cuda-backend
# Plus a final cluster smoke test:
bash deploy_cuda.sh
ssh mimo "tail -200 /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda/cuda_test_*.out | grep -E 'PASS|FAIL'"
```

- [ ] **Step 6: Commit the checklist**

```bash
git add docs/superpowers/reviews/2026-04-27-cuda-backend-merge-checklist.md
git commit -m "Add cuda-backend → master merge checklist with rebase plan"
```

- [ ] **Step 7: Hand back to user for the actual rebase**

The rebase is destructive history rewriting on a branch the user has been pushing — do NOT execute it. Print the recommended `git rebase -i` invocation and let the user run it.

---

## Self-Review

- **Spec coverage:** All 5 punch-list items have a task. ✓
- **Placeholder scan:** No "TBD"/"implement later"; every command is a real command. ✓
- **Type/identifier consistency:** Workspace member naming follows the Perf #1 precedent (`*_workspace`). Subagent names match `.claude/agents/` files. Test names match `test/cuda/`. ✓
