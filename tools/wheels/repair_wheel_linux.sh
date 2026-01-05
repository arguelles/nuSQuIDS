#!/bin/bash
# Repair Linux wheel: fix RPATH and then run auditwheel
# Usage: repair_wheel_linux.sh <wheel> <dest_dir>

set -ex

WHEEL="$1"
DEST_DIR="$2"

# Create temp directory
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Unpack the wheel
unzip -q "$WHEEL" -d "$TMPDIR"

# Find all .so files in the nuSQuIDS directory
PACKAGE_DIR="$TMPDIR/nuSQuIDS"

# Fix RPATH for all shared libraries
fix_rpath() {
    local BINARY="$1"
    if [ ! -f "$BINARY" ]; then
        return
    fi

    echo "Setting RPATH in: $BINARY"
    patchelf --set-rpath '$ORIGIN' "$BINARY" 2>/dev/null || true
}

# Fix all .so files in the package
for f in "$PACKAGE_DIR"/*.so "$PACKAGE_DIR"/*.so.*; do
    if [ -f "$f" ]; then
        fix_rpath "$f"
    fi
done

# Also check the lib directory if it exists
if [ -d "$TMPDIR/lib" ]; then
    for f in "$TMPDIR/lib"/*.so "$TMPDIR/lib"/*.so.*; do
        if [ -f "$f" ]; then
            fix_rpath "$f"
        fi
    done
fi

# Repack the wheel
WHEEL_NAME=$(basename "$WHEEL")
cd "$TMPDIR"
zip -q -r "$TMPDIR/fixed_$WHEEL_NAME" ./*

# Run auditwheel on the fixed wheel
auditwheel repair -w "$DEST_DIR" "$TMPDIR/fixed_$WHEEL_NAME"
