#!/bin/bash
# Repair macOS wheel: fix library paths and then run delocate
# Usage: repair_wheel_macos.sh <wheel> <dest_dir> <delocate_archs>

set -ex

WHEEL="$1"
DEST_DIR="$2"
DELOCATE_ARCHS="$3"

# Create temp directory
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Unpack the wheel
unzip -q "$WHEEL" -d "$TMPDIR"

# Find all dylibs and .so files in the nuSQuIDS directory
PACKAGE_DIR="$TMPDIR/nuSQuIDS"

fix_lib_refs() {
    local BINARY="$1"
    if [ ! -f "$BINARY" ]; then
        return
    fi

    echo "Fixing library paths in: $BINARY"

    # Get current dependencies
    DEPS=$(otool -L "$BINARY" | tail -n +2)

    # Fix all absolute paths to use @loader_path
    echo "$DEPS" | while read -r line; do
        # Extract the path (first field before the parenthesis)
        DEP_PATH=$(echo "$line" | sed -E 's/^[[:space:]]*([^ ]+).*/\1/')

        # Check if it's an absolute path to a library we care about
        if [[ "$DEP_PATH" == /*/libSQuIDS* ]] || [[ "$DEP_PATH" == /*/libnuSQuIDS* ]]; then
            LIB_NAME=$(basename "$DEP_PATH")
            echo "  Changing $DEP_PATH -> @loader_path/$LIB_NAME"
            install_name_tool -change "$DEP_PATH" "@loader_path/$LIB_NAME" "$BINARY" 2>/dev/null || true
        fi
    done

    # Also fix common non-absolute references
    install_name_tool -change "libSQuIDS.1.dylib" "@loader_path/libSQuIDS.1.dylib" "$BINARY" 2>/dev/null || true
    install_name_tool -change "libnuSQuIDS.1.dylib" "@loader_path/libnuSQuIDS.1.dylib" "$BINARY" 2>/dev/null || true
    install_name_tool -change "lib/nuSQuIDS.so" "@loader_path/nuSQuIDS.so" "$BINARY" 2>/dev/null || true
}

# Fix all binaries in the package
for f in "$PACKAGE_DIR"/*.so "$PACKAGE_DIR"/*.dylib; do
    if [ -f "$f" ]; then
        fix_lib_refs "$f"

        # Set the install name ID for dylibs
        if [[ "$f" == *.dylib ]]; then
            LIB_NAME=$(basename "$f")
            install_name_tool -id "@loader_path/$LIB_NAME" "$f" 2>/dev/null || true
        fi
    fi
done

# Also fix binaries in the lib directory if it exists
if [ -d "$TMPDIR/lib" ]; then
    for f in "$TMPDIR/lib"/*.so "$TMPDIR/lib"/*.dylib; do
        if [ -f "$f" ]; then
            fix_lib_refs "$f"
            # Set the install name ID for dylibs
            if [[ "$f" == *.dylib ]]; then
                LIB_NAME=$(basename "$f")
                install_name_tool -id "@loader_path/$LIB_NAME" "$f" 2>/dev/null || true
            fi
        fi
    done
fi

# Repack the wheel with original name (must match dist-info directory)
WHEEL_NAME=$(basename "$WHEEL")
cd "$TMPDIR"
zip -q -r "$TMPDIR/$WHEEL_NAME" ./*

# Run delocate on the fixed wheel
# Use --ignore-missing-dependencies to handle libraries that may not be found in expected locations
# Use --check-archs warn to allow version mismatches (Homebrew libs may target newer macOS)
delocate-wheel --require-archs "$DELOCATE_ARCHS" -w "$DEST_DIR" -v --ignore-missing-dependencies "$TMPDIR/$WHEEL_NAME"
