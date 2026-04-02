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

source ../env_vars.sh

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
echo "=== Compiling benchmark_cuda ==="
$CXX $CXXFLAGS $CFLAGS -o benchmark_cuda benchmark_cuda.cpp $LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1
echo "=== Running benchmark_cuda ==="
./benchmark_cuda
