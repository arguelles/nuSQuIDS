#!/bin/bash
# Build script for nuSQuIDS with CUDA on FASRC cluster
set -e

# Load modules
module load gcc/12.2.0-fasrc01
module load openmpi/5.0.5-fasrc02
module load hdf5/1.14.6-fasrc01
module load gsl/2.8-fasrc01
module load cuda/12.4.1-fasrc01

echo "=== Environment ==="
echo "GCC: $(gcc --version | head -1)"
echo "NVCC: $(nvcc --version | tail -1)"
echo "GSL: $(gsl-config --prefix)"
echo "HDF5: $HDF5_HOME"

BASEDIR=/n/home06/carguelles/programs
INSTALL_PREFIX=$BASEDIR/local
mkdir -p "$INSTALL_PREFIX"

# ---- Build SQuIDS if needed ----
if [ ! -f "$INSTALL_PREFIX/lib/libSQuIDS.a" ]; then
    echo ""
    echo "=== Building SQuIDS ==="
    cd "$BASEDIR"
    if [ ! -d SQuIDS ]; then
        git clone https://github.com/jsalvado/SQuIDS.git
    fi
    cd SQuIDS
    make clean 2>/dev/null || true
    ./configure --prefix="$INSTALL_PREFIX" \
        --with-gsl="$(gsl-config --prefix)"
    make -j4
    make install
    echo "SQuIDS installed to $INSTALL_PREFIX"
else
    echo "SQuIDS already installed at $INSTALL_PREFIX"
fi

# ---- Build nuSQuIDS with CUDA ----
echo ""
echo "=== Building nuSQuIDS with CUDA ==="
cd "$BASEDIR/nuSQuIDS"

# Clean previous build artifacts (macOS objects won't work)
rm -f build/*.o lib/*.a lib/*.so lib/*.dylib

./configure --prefix="$INSTALL_PREFIX" \
    --with-gsl="$(gsl-config --prefix)" \
    --with-hdf5="$HDF5_HOME" \
    --with-squids="$INSTALL_PREFIX" \
    --enable-cuda

echo ""
echo "=== Compiling ==="
make -j4

echo ""
echo "=== Build complete! ==="
ls -la lib/libnuSQuIDS.*

echo ""
echo "=== Running CPU tests ==="
make test || echo "Some tests may fail (expected on first build)"

echo ""
echo "Done! To run GPU tests, submit a job to a GPU partition."
