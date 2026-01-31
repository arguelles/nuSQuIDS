#!/usr/bin/env python3
"""
Generate SHA256 hashes for nuSQuIDS data files.

This script generates the hash registry needed for pooch to verify data integrity.
Run this whenever data files are updated and paste the output into nuSQuIDS_data.py.

Usage:
    python tools/generate_data_hashes.py
"""

import hashlib
from pathlib import Path


def sha256_hash(filepath):
    """Calculate SHA256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def main():
    # Find the data directory
    script_dir = Path(__file__).parent
    data_dir = script_dir.parent / "data"

    if not data_dir.exists():
        print(f"Error: Data directory not found at {data_dir}")
        return

    print("# Registry of data files with SHA256 hashes")
    print("# Copy this into nuSQuIDS_data.py REGISTRY dict")
    print("REGISTRY = {")

    # Process xsections
    xsections_dir = data_dir / "xsections"
    if xsections_dir.exists():
        print("    # Cross-section files")
        for ext in ["*.h5", "*.dat"]:
            for f in sorted(xsections_dir.glob(ext)):
                if f.name == "README.md":
                    continue
                rel_path = f"xsections/{f.name}"
                file_hash = sha256_hash(f)
                print(f'    "{rel_path}": "{file_hash}",')

    # Process atmos_prod
    atmos_dir = data_dir / "atmos_prod"
    if atmos_dir.exists():
        print("    # Atmospheric production profiles")
        for f in sorted(atmos_dir.glob("*.dat")):
            rel_path = f"atmos_prod/{f.name}"
            file_hash = sha256_hash(f)
            print(f'    "{rel_path}": "{file_hash}",')

    print("}")

    # Also print total size
    print("\n# File sizes:")
    total_size = 0
    for subdir in ["xsections", "atmos_prod"]:
        subdir_path = data_dir / subdir
        if subdir_path.exists():
            for f in subdir_path.iterdir():
                if f.is_file() and f.name != "README.md":
                    size = f.stat().st_size
                    total_size += size
                    print(f"# {subdir}/{f.name}: {size / 1024 / 1024:.1f} MB")

    print(f"# Total: {total_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
