# ODE solver benchmarks

Tools to compare the integration methods available to nuSQuIDS: the GSL
steppers (the historical path), the native embedded Runge-Kutta pairs of SQuIDS,
per-node integration, and the exponential (Magnus-type) node integrators.

Methods are compared at **equal achieved accuracy**, not at equal nominal
tolerance: every run is scored by its maximum absolute error in the flavor
expectation values against a converged reference, and for a target accuracy
the cost of a method is its cheapest run that reached it.

## Building and running

```bash
./build.sh                                   # uses the flags ./configure wrote into ../../Makefile
./ode_bench <workload> <method> <tol> [out.txt] [reps]
```

`ode_bench` prints the time spent in `EvolveState` (setup excluded), the number
of right-hand-side evaluations and the accepted/rejected steps. For per-node
methods the evaluation count is converted to full-system equivalents.

References (written to `ref/`, not in git):

```bash
for w in earth_osc const_osc atm_osc earth_int sterile earth_subgev earth_mev; do
  ./ode_bench $w rk8pd 1e-13 ref/${w}_rk8pd_1e-13.txt
  ./ode_bench $w rk8pd 1e-14 ref/${w}_rk8pd_1e-14.txt   # the reference
done
```

Tightening rk8pd from 1e-13 to 1e-14 moves every workload by at most 6e-12,
and for constant density the reference agrees with the exact solution to 4e-14,
so errors down to ~1e-11 are resolved.

Work-precision sweep and analysis:

```bash
./snapshot.sh /tmp/frozen                    # freeze binary + libraries for long sweeps
BENCH=/tmp/frozen/ode_bench ./sweep.sh out.csv 3 "earth_osc atm_osc" "rkf45 n_vern7:I" "1e-4 1e-6 1e-8"
python3 frontier.py out.csv                  # fastest methods per target accuracy
python3 rank.py out.csv                      # geometric-mean speedup at equal accuracy
python plot_wp.py out.csv wp.pdf             # work-precision diagrams (matplotlib)
```

## Workloads

| name | problem |
|---|---|
| `earth_osc` | core-crossing chord (12742 km, PREM), 0.5-500 GeV, 150 nodes, nu and nubar |
| `const_osc` | same energies through constant density (exact solution exists) |
| `atm_osc` | `nuSQUIDSAtm`, 20 zeniths, 1-1000 GeV, 100 nodes |
| `earth_int` | core-crossing, 0.1-10^4 TeV, NC interactions and tau regeneration |
| `sterile` | 3+1, dm2_41 = 1 eV^2, theta_24 = 0.1, 0.1-100 TeV (Earth MSW resonance) |
| `earth_subgev` | core-crossing, 0.1-2 GeV |
| `earth_mev` | core-crossing, 2-20 MeV (solar day-night) |

Oscillation parameters: NuFIT 5.2 (2022), normal ordering, without SK
atmospheric data.

## Method names

| name | meaning |
|---|---|
| `rkf45`, `rkck`, `rk8pd`, `rk4`, `msadams` | GSL steppers through the GSL driver (`rkf45` is the default) |
| `n_<m>` | native SQuIDS method `<m>` (`rkf45 bs3 dp5 tsit5 vern6 vern7 vern8 vern9 rkv76 dop853`), PI controller |
| `n_<m>:I`, `n_<m>:H211b` | the same with the I or H211b controller |
| `ni_<m>[:ctl]` | each energy node integrated separately with native method `<m>` |
| `ni_magnus4`, `ni_cfmagnus4`, `ni_magnus6`, `ni_modmagnus` | exponential node integrators |
| `analytic` | nuSQuIDS' exact constant-density evolution |

## Results

Measured on an Apple M-series CPU (single thread), 3 repetitions per point,
7 tolerances from 1e-4 to 1e-10. Speedups are geometric means over target
accuracies 1e-5, 1e-6, 1e-7 and 1e-8, at equal achieved accuracy.

**Earth propagation, GeV to PeV** (`earth_osc`, `atm_osc`, `earth_int`, `sterile`):

| method | vs GSL rkf45 (default) | vs best GSL stepper |
|---|---|---|
| native Vern9, I controller | 2.08x | 1.64x |
| native Vern7, I controller | 2.04x | 1.62x |
| native DOP853, I | 1.96x | 1.55x |
| native DP5, I | 1.79x | 1.42x |
| native RKF45, I (same method as the default) | 1.56x | 1.23x |

The last row isolates the driver: the default's own Fehlberg coefficients run
through the native driver (first-same-as-last, warm start, unit-free initial
step, stops at PREM breakpoints) gain 1.56x; the higher-order method adds the
rest. Per workload, Vern7/Vern9 are 1.5-2.4x faster on the core-crossing chord,
1.3-3.4x on the atmospheric grid and 2-7x with interactions. On the atmospheric
grid the default did not reach 1e-8 at any tolerance down to 1e-10.

**Low energies** (`earth_subgev`, `earth_mev`): native Vern9 is 1.8x faster
than the default. At MeV, the modified Magnus node integrator takes ~300 km
steps where Runge-Kutta needs ~2 km and is 8.6-17x faster down to 1e-6
accuracy; it does not reach 1e-7 (the neglected second Magnus term), where
Vern9 is the fastest (2-2.5x).

**Constant density**: the Magnus node integrators are exact (Magnus6 20x,
Magnus4 14x faster than the default at equal accuracy), and per-node Vern7 is
5x faster; nuSQuIDS' analytic shortcut (`Set_AllowConstantDensityOscillationOnlyEvolution`)
remains the best choice there.

What did not pay off: the PI and H211b controllers (the I controller needed
less work on every Earth workload); GSL's multistep `msadams` (fewer
evaluations but no less time); per-node Runge-Kutta on Earth profiles (every
node meets every density feature); standard Magnus methods above ~100 MeV
(their error grows with |H0|, so they take as many steps as Runge-Kutta at a
much higher cost per step). Implicit methods were not tried: SQuIDS provides
no Jacobian, and these problems are not stiff (absorption rate times step size
stays far inside the explicit stability region; rejections stay low at PeV).
Extrapolation and Taylor methods were not tried either: the first needs a
smooth right-hand side (the Akima spline of PREM is only C1), the second
derivatives of the whole right-hand side.

Recommended settings:

```cpp
nus.Set_ODEMethod(squids::ODEMethod::Vern7);              // GeV and above
nus.Set_NodeEvolution(nusquids::NodeEvolution::ModifiedMagnus); // MeV, no interactions
```

