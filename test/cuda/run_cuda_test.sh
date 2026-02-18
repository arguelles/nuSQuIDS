#!/bin/bash
#SBATCH -J nusquids_cuda_test
#SBATCH -p arguelles_delgado_gpu_a100
#SBATCH --gres=gpu:1
#SBATCH -t 00:30:00
#SBATCH --mem=8G
#SBATCH -o cuda_test_%j.out
#SBATCH -e cuda_test_%j.err

# Use gsl/2.7 to match libSQuIDS.so (built against libgsl.so.25)
module load gcc/12.2.0 cuda/12.4.1 openmpi/5.0.5 hdf5/1.14.6 gsl/2.7

cd /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda

source ../env_vars.sh

# Override library paths to use correct GSL 2.7 and HDF5 1.14.6
GSL27_LIBDIR="/n/sw/helmod-rocky8/apps/Core/gsl/2.7-fasrc01/lib64"
HDF5_LIBDIR="/n/sw/helmod-rocky8/apps/MPI/gcc/12.2.0-fasrc01/openmpi/5.0.5-fasrc02/hdf5/1.14.6-fasrc01/lib"
export LD_LIBRARY_PATH="${GSL27_LIBDIR}:${HDF5_LIBDIR}:${LD_LIBRARY_PATH}"
EXTRA_LDFLAGS="-L${GSL27_LIBDIR} -Wl,-rpath,${GSL27_LIBDIR} -L${HDF5_LIBDIR} -Wl,-rpath,${HDF5_LIBDIR}"

echo "GSL 2.7 lib: ${GSL27_LIBDIR}"
ls -la ${GSL27_LIBDIR}/libgsl.so.25 2>&1
echo ""

echo "=== Compiling cuda_test ==="
$CXX $CXXFLAGS $CFLAGS -o cuda_test cuda_test.cpp $LDFLAGS $EXTRA_LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1
echo "=== Running cuda_test ==="
./cuda_test

echo ""
echo "=== Compiling benchmark_cuda ==="
$CXX $CXXFLAGS $CFLAGS -o benchmark_cuda benchmark_cuda.cpp $LDFLAGS $EXTRA_LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1
echo "=== Running benchmark_cuda ==="
./benchmark_cuda
