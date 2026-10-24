#!/bin/bash
# Freeze the benchmark binary together with the exact libnuSQuIDS and libSQuIDS
# it was built against, so long sweeps are not disturbed by later rebuilds.
#
# usage: snapshot.sh <dir>     (then run <dir>/ode_bench, or sweep.sh with BENCH=<dir>/ode_bench)
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
dest="$1"
mkdir -p "$dest"
cp "$here/ode_bench" "$dest/ode_bench"
cp "$root/lib/libnuSQuIDS.1.dylib" "$dest/libnuSQuIDS.1.dylib"
squids_lib="$(otool -L "$root/lib/libnuSQuIDS.1.dylib" | awk '/libSQuIDS/{print $1}')"
cp "$squids_lib" "$dest/libSQuIDS.1.dylib"
# short relative names always fit in the existing load commands
install_name_tool -change "$root/lib/libnuSQuIDS.1.dylib" "@executable_path/libnuSQuIDS.1.dylib" "$dest/ode_bench"
install_name_tool -change "$squids_lib" "@executable_path/libSQuIDS.1.dylib" "$dest/ode_bench"
install_name_tool -change "$squids_lib" "@loader_path/libSQuIDS.1.dylib" "$dest/libnuSQuIDS.1.dylib"
codesign -f -s - "$dest/ode_bench" "$dest/libnuSQuIDS.1.dylib" "$dest/libSQuIDS.1.dylib" 2>/dev/null || true
otool -L "$dest/ode_bench" | grep -i squids
echo "snapshot in $dest"
