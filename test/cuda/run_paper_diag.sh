#!/bin/bash
# Run the paper-diagnostics binary. Defaults target the lab MIG partition;
# override with `sbatch -p gpu --gres=gpu:nvidia_a100-sxm4-80gb:1 ...` to
# pin a full A100 from the public `gpu` partition.
#SBATCH -J nusquids_paper_diag
#SBATCH -p arguelles_delgado_gpu_a100
#SBATCH --gres=gpu:1
#SBATCH -t 01:00:00
#SBATCH --mem=16G
#SBATCH -o paper_diag_%j.out
#SBATCH -e paper_diag_%j.err

module load gcc/12.2.0 cuda/12.4.1 openmpi/5.0.5 hdf5/1.14.6 gsl/2.8

cd /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/nuSQuIDS/test/cuda

source ../env_vars.sh

echo "=== GPU info ==="
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv

echo ""
echo "=== Compiling cuda_paper_diagnostics ==="
$CXX $CXXFLAGS $CFLAGS -o cuda_paper_diagnostics cuda_paper_diagnostics.cpp $LDFLAGS -L$CUDA_HOME/lib64/stubs -lcuda 2>&1
if [ $? -ne 0 ]; then
  echo "=== COMPILATION FAILED ==="
  exit 1
fi

echo ""
echo "=== Running cuda_paper_diagnostics ==="
./cuda_paper_diagnostics

echo ""
echo "=== Output files ==="
ls -la paper_diag_*${SLURM_JOB_ID}*
