#!/bin/bash
# Deploy and test CUDA backend on FASRC cluster
# Run this locally: bash deploy_cuda.sh

set -e

REMOTE_HOST="mimo"
REMOTE_DIR="/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS"

echo "=== Step 1: Syncing code to cluster ==="
rsync -avz --delete \
  --exclude='build/' \
  --exclude='lib/' \
  --exclude='test/env_vars.sh' \
  --exclude='.git/' \
  --exclude='*.o' \
  --exclude='*.so' \
  --exclude='*.dylib' \
  --exclude='*.a' \
  --exclude='deploy_cuda.sh' \
  /Users/carguelles/Programs/nuSQuIDS/ \
  ${REMOTE_HOST}:${REMOTE_DIR}/

echo ""
echo "=== Step 2: Building on cluster ==="
ssh ${REMOTE_HOST} "bash -l -c '
  cd ${REMOTE_DIR}
  module load gcc/12.2.0 cuda/12.4.1 openmpi/5.0.5 hdf5/1.14.6 gsl/2.8

  SQUIDS_PATH=/n/home06/carguelles/programs/local
  echo \"Using SQuIDS at: \${SQUIDS_PATH}\"

  echo \"\"
  echo \"--- Configuring ---\"
  CXX=g++ ./configure \
    --with-squids=\${SQUIDS_PATH} \
    --with-hdf5=\${HDF5_HOME} \
    --enable-cuda

  echo \"\"
  echo \"--- Building (clean) ---\"
  make clean
  mkdir -p build lib
  make -j4

  echo \"\"
  echo \"--- Build complete ---\"
  ls -la lib/
'"

echo ""
echo "=== Step 3: Submitting test job ==="
ssh ${REMOTE_HOST} "bash -l -c '
  cd ${REMOTE_DIR}/test/cuda
  sbatch run_cuda_test.sh
'"

echo ""
echo "=== Done! ==="
echo "Monitor job with: ssh ${REMOTE_HOST} 'squeue -u \$USER'"
echo "Check output with: ssh ${REMOTE_HOST} 'cat ${REMOTE_DIR}/test/cuda_test_*.out'"
