#!/bin/bash
# Script to install dependencies before wheel build
# Called by cibuildwheel via before-all hook

set -xe

PROJECT_DIR="$1"

echo "Project directory: $PROJECT_DIR"
echo "Runner OS: $RUNNER_OS"

if [[ $RUNNER_OS == "Linux" ]]; then
    # manylinux_2_28 uses AlmaLinux 8 with dnf
    if command -v dnf &> /dev/null; then
        # Enable EPEL for hdf5-devel and patchelf
        dnf install -y epel-release
        dnf install -y hdf5-devel gsl-devel patchelf zip
    elif command -v yum &> /dev/null; then
        # For older manylinux, enable EPEL
        yum install -y epel-release
        yum install -y hdf5-devel gsl-devel patchelf zip
    elif command -v apk &> /dev/null; then
        # musllinux uses Alpine
        apk add hdf5-dev gsl-dev patchelf zip
    else
        echo "Unknown Linux package manager" 1>&2
        exit 1
    fi
elif [[ $RUNNER_OS == "macOS" ]]; then
    brew install hdf5 gsl
else
    echo "Unknown runner OS: $RUNNER_OS" 1>&2
    exit 1
fi

echo "Dependencies installed successfully"
