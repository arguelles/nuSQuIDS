#!/bin/bash
# Build the ODE benchmark against this checkout, using the compiler flags that
# ./configure wrote into the top-level Makefile.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
cd "$root"
vars="$(mktemp)"
trap 'rm -f "$vars"' EXIT
printf 'include Makefile\nprint-%%:\n\t@echo $($*)\n' > "$vars"
CXX="$(make -s -f "$vars" print-CXX)"
FLAGS="$(make -s -f "$vars" print-EXAMPLES_FLAGS)"
LDFLAGS="$(make -s -f "$vars" print-LDFLAGS)"
# shellcheck disable=SC2086
pad=""
[ "$(uname -s)" = Darwin ] && pad="-Wl,-headerpad_max_install_names"
$CXX $FLAGS ${EXTRA_FLAGS:-} "$here/ode_bench.cpp" -lnuSQuIDS $LDFLAGS $pad -o "$here/ode_bench"
# libnuSQuIDS has a bare install name; pin it so the binary runs without
# DYLD_LIBRARY_PATH (which macOS strips from child processes of /bin/bash)
if [ "$(uname -s)" = Darwin ]; then
  install_name_tool -change libnuSQuIDS.1.dylib "$root/lib/libnuSQuIDS.1.dylib" "$here/ode_bench"
fi
echo "built $here/ode_bench"
