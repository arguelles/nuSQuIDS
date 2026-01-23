#!/bin/bash
# Build script for nuSQuIDS Julia bindings
#
# Prerequisites:
#   1. Julia must be installed with CxxWrap.jl package
#   2. nuSQuIDS must be installed (headers in /usr/local/include, library in /usr/local/lib)
#   3. libcxxwrap-julia must be available
#
# Usage:
#   ./build.sh
#
# The script will:
#   1. Get the libcxxwrap-julia prefix from Julia
#   2. Build the C++ wrapper library using CMake
#   3. Place the library in resources/julia/lib/

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
LIB_DIR="${SCRIPT_DIR}/lib"

echo "=== Building nuSQuIDS Julia Bindings ==="

# Check for Julia
if ! command -v julia &> /dev/null; then
    echo "Error: Julia not found. Please install Julia first."
    exit 1
fi

# Check that CxxWrap is installed
echo "Checking for CxxWrap.jl..."
CXXWRAP_PREFIX=$(julia -e 'using CxxWrap; print(CxxWrap.prefix_path())' 2>/dev/null)
if [ -z "$CXXWRAP_PREFIX" ]; then
    echo "Error: CxxWrap.jl not found. Please install it:"
    echo "  julia -e 'using Pkg; Pkg.add(\"CxxWrap\")'"
    exit 1
fi
echo "Found CxxWrap at: $CXXWRAP_PREFIX"

# Get Julia include directory
JULIA_DIR=$(julia -e 'print(dirname(Sys.BINDIR))')
echo "Julia installation: $JULIA_DIR"

# Create build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure with CMake
echo "Configuring CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$CXXWRAP_PREFIX" \
    -DJulia_PREFIX="$JULIA_DIR"

# Build
echo "Building..."
cmake --build . --config Release -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)

# Create lib directory and copy library
mkdir -p "$LIB_DIR"
if [ -f "$BUILD_DIR/nuSQuIDSjl.dylib" ]; then
    cp "$BUILD_DIR/nuSQuIDSjl.dylib" "$LIB_DIR/"
elif [ -f "$BUILD_DIR/nuSQuIDSjl.so" ]; then
    cp "$BUILD_DIR/nuSQuIDSjl.so" "$LIB_DIR/"
elif [ -f "$LIB_DIR/nuSQuIDSjl.dylib" ] || [ -f "$LIB_DIR/nuSQuIDSjl.so" ]; then
    echo "Library already in lib directory"
else
    echo "Warning: Could not find built library"
fi

echo ""
echo "=== Build Complete ==="
echo ""
echo "To use the Julia bindings:"
echo "  1. Add the package to your Julia environment:"
echo "     julia> ] develop ${SCRIPT_DIR}/nuSQuIDS"
echo ""
echo "  2. Then use it:"
echo "     julia> using nuSQuIDS"
echo ""
