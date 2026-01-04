#!/bin/bash
# Script to install dependencies before wheel build
# Called by cibuildwheel via before-all hook

set -xe

PROJECT_DIR="$1"

echo "Project directory: $PROJECT_DIR"
echo "Runner OS: $RUNNER_OS"

if [[ $RUNNER_OS == "Linux" ]]; then
    # manylinux_2_28 uses AlmaLinux 8 with dnf/yum
    if command -v dnf &> /dev/null; then
        dnf install -y hdf5-devel gsl-devel
    elif command -v yum &> /dev/null; then
        yum install -y hdf5-devel gsl-devel
    elif command -v apk &> /dev/null; then
        # musllinux uses Alpine
        apk add hdf5-dev gsl-dev
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
