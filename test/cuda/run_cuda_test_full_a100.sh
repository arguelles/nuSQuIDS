#!/bin/bash
# Full A100-SXM4-80GB run — uses the public `gpu` partition since the
# lab partition (arguelles_delgado_gpu_a100) is MIG-only. Keeps the same
# build/run flow as run_cuda_test.sh.
#SBATCH -J nusquids_cuda_test_full_a100
#SBATCH -p gpu
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH -t 00:30:00
#SBATCH --mem=16G
#SBATCH -o cuda_test_full_a100_%j.out
#SBATCH -e cuda_test_full_a100_%j.err

module load gcc/12.2.0 cuda/12.4.1 openmpi/5.0.5 hdf5/1.14.6 gsl/2.8

cd /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda

source ../env_vars.sh

echo "=== GPU info ==="
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv

echo ""
echo "=== Compiling cuda_rhs_diagnostic ==="
$CXX $CXXFLAGS $CFLAGS -o cuda_rhs_diagnostic cuda_rhs_diagnostic.cpp $LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1
if [ $? -eq 0 ]; then
  echo "=== Running cuda_rhs_diagnostic ==="
  ./cuda_rhs_diagnostic
else
  echo "=== DIAGNOSTIC COMPILATION FAILED ==="
fi

echo ""
echo "=== Compiling cuda_test ==="
$CXX $CXXFLAGS $CFLAGS -o cuda_test cuda_test.cpp $LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1
echo "=== Running cuda_test ==="
./cuda_test

echo ""
echo "=== Compiling cuda_single_baseline_test ==="
$CXX $CXXFLAGS $CFLAGS -o cuda_single_baseline_test cuda_single_baseline_test.cpp $LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1
echo "=== Running cuda_single_baseline_test ==="
./cuda_single_baseline_test

echo ""
echo "=== Compiling cuda_interactions_test ==="
$CXX $CXXFLAGS $CFLAGS -o cuda_interactions_test cuda_interactions_test.cpp $LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1
if [ $? -eq 0 ]; then
  echo "=== Running cuda_interactions_test ==="
  ./cuda_interactions_test
else
  echo "=== COMPILATION FAILED ==="
fi

echo ""
echo "=== Compiling cuda_dopri5_fsal_test ==="
$CXX $CXXFLAGS $CFLAGS -o cuda_dopri5_fsal_test cuda_dopri5_fsal_test.cpp $LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1
if [ $? -eq 0 ]; then
  echo "=== Running cuda_dopri5_fsal_test ==="
  ./cuda_dopri5_fsal_test
else
  echo "=== DOPRI5 FSAL TEST COMPILATION FAILED ==="
fi

echo ""
echo "=== Compiling benchmark_cuda ==="
$CXX $CXXFLAGS $CFLAGS -o benchmark_cuda benchmark_cuda.cpp $LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1
echo "=== Running benchmark_cuda ==="
./benchmark_cuda
