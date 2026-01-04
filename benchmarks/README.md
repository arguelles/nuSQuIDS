# nuSQuIDS Benchmark Suite

This directory contains performance benchmarks for nuSQuIDS.

## Running Benchmarks

From the project root:

```bash
# Full benchmark (more iterations, more accurate)
make benchmark

# Quick benchmark (fewer iterations, faster)
make benchmark-quick
```

## Benchmark Tests

The benchmark suite measures performance across different propagation modes and physics configurations:

| Test | Description |
|------|-------------|
| Single Energy (isoscalar) | Single nu_mu through Earth (100 GeV) |
| Multi-E isoscalar (no int.) | Power-law spectrum, 1 GeV - 10 TeV, 200 nodes |
| Multi-E isoscalar (with int.) | Power-law spectrum, 10 GeV - 1 PeV, 100 nodes, NC+CC |
| Atm isoscalar (no int.) | nuSQUIDSAtm, 50 energies x 20 zenith angles |
| Atm isoscalar (with int.) | nuSQUIDSAtm with NC+CC interactions |
| Atm isoscalar (int.+Glashow) | Adds Glashow resonance (extends to 10 PeV) |
| Atm isoscalar (full physics) | NC+CC+Glashow+Tau regeneration |
| Atm nuclear (PREM composition) | Uses PREM with element composition |

## Example Output

```
================================================================================
                        nuSQuIDS Performance Benchmark
                              Version: 1.13.0
================================================================================

System Info:
  Platform: macOS (arm64)
  Date: Sat Jan  3 12:00:00 2026

--------------------------------------------------------------------------------
Single Energy Mode (isoscalar)
--------------------------------------------------------------------------------
  Test: Propagate nu_mu through Earth (cos_zen = -1), E = 100 GeV
  Iterations: 100

  Time per iteration:  4.52 ms
  Total time:          452 ms
  Status: PASS

...

================================================================================
                               Summary
================================================================================
  Single Energy (isoscalar)        4.52 ms  [PASS]
  Multi-E isoscalar (no int.)      85.3 ms  [PASS]
  Multi-E isoscalar (with int.)    1.25 s   [PASS]
  ...

  Total benchmark time: 45.2 s
  Tests passed: 8/8
================================================================================
```

## Command Line Options

```bash
# Quick mode (fewer iterations)
./benchmark --quick
./benchmark -q
```

## Using Benchmarks for Performance Regression

Run benchmarks before and after code changes to detect performance regressions:

```bash
# Before changes
make benchmark > baseline.txt

# After changes
make benchmark > current.txt

# Compare results
diff baseline.txt current.txt
```

## Adding New Benchmarks

To add a new benchmark, create a function following this pattern:

```cpp
BenchmarkResult benchmark_my_test(int iterations = 10) {
    BenchmarkResult result;
    result.name = "My Test Name";
    result.description = "Description of what is being tested";
    result.iterations = iterations;
    result.passed = true;

    auto start = Clock::now();

    for (int iter = 0; iter < iterations; iter++) {
        // ... benchmark code ...
    }

    auto end = Clock::now();
    result.total_time_ms = Duration(end - start).count();
    result.time_per_iter_ms = result.total_time_ms / iterations;

    return result;
}
```

Then add it to `main()`:

```cpp
print_test_header("My Test Category");
auto r = benchmark_my_test(quick_mode ? 1 : 10);
print_test_result(r);
results.push_back(r);
```
