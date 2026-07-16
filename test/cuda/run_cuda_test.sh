#!/bin/bash
#SBATCH -J nusquids_cuda_test
#SBATCH -p arguelles_delgado_gpu_a100
#SBATCH --gres=gpu:1
#SBATCH -t 00:30:00
#SBATCH --mem=8G
#SBATCH -o cuda_test_%j.out
#SBATCH -e cuda_test_%j.err

module load gcc/12.2.0 cuda/12.4.1 openmpi/5.0.5 hdf5/1.14.6 gsl/2.8

cd /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda

# Run the standard CUDA test suite (compile + run each test, produce Report.txt).
# Fails non-zero if any test fails; the SLURM job status will reflect this.
echo "=== Running cuda test suite via run_cuda_tests ==="
bash run_cuda_tests
suite_rc=$?

# The benchmark and rhs diagnostic are not part of the standard suite; they
# produce performance / debug output rather than PASS/FAIL. Run them
# separately so their logs continue to land in cuda_test_${SLURM_JOB_ID}.out.
source ../env_vars.sh
echo ""
echo "=== Compiling and running cuda_rhs_diagnostic (diagnostic; not gated) ==="
if $CXX $CXXFLAGS $CFLAGS -o cuda_rhs_diagnostic cuda_rhs_diagnostic.cpp \
    $LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1; then
  ./cuda_rhs_diagnostic
else
  echo "=== DIAGNOSTIC COMPILATION FAILED ==="
fi

echo ""
echo "=== Compiling and running benchmark_cuda (perf sweep; not gated) ==="
if $CXX $CXXFLAGS $CFLAGS -o benchmark_cuda benchmark_cuda.cpp \
    $LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1; then
  ./benchmark_cuda
else
  echo "=== BENCHMARK COMPILATION FAILED ==="
fi

# Propagate the suite exit code to SLURM.
exit $suite_rc
