"""
nuSQuIDS Data Manager

This module handles downloading and caching of large data files (cross sections,
atmospheric production profiles) that cannot be included in the PyPI package due
to size limits.

Data files are downloaded from GitHub releases and cached in the user's home directory.
"""

import os
import sys
import hashlib
from pathlib import Path

try:
    import pooch
    HAS_POOCH = True
except ImportError:
    HAS_POOCH = False

# Version of the data files (update when data files change)
DATA_VERSION = "v1.13.3"

# Base URL for data files (GitHub release assets)
GITHUB_REPO = "arguelles/nuSQuIDS"
BASE_URL = f"https://github.com/{GITHUB_REPO}/releases/download/{DATA_VERSION}/"

# Alternative: use a dedicated data repository
# BASE_URL = "https://github.com/arguelles/nuSQuIDS-data/raw/main/"

# Registry of data files with their SHA256 hashes
# Generate hashes with: python tools/generate_data_hashes.py
# These should be updated when data files change
REGISTRY = {
    # CSMS cross sections (~51 MB total)
    "xsections/csms.h5": "1a475aef6dee2f8c777cc9568407231bdaa0dbf0224bc33fbd880c55cbf9daa1",
    "xsections/csms_proton.h5": "38f952442061849a5544f47063e8bb590b692f9343a7b57c69e938401e358d4c",
    "xsections/csms_neutron.h5": "d03eb47a9ddcb41aa77eb727f90e848e23415380b5871cb192a0c6b5188fa210",
    "xsections/csms_square.h5": "e749c6f39d1fa3bf42247be45a9ff731a241a01a2fc8edabc211fb4876f69d57",

    # NuSigma cross sections (~30 MB total)
    "xsections/nusigma_sigma_CC.dat": "81e90d5e8bf98eaac309c61e61d5a4f1b6b00b9ba479b42990106beaf186598f",
    "xsections/nusigma_sigma_NC.dat": "543b2f31b33fc8b4b8026475fe2cc07f04daf02e207cabe993c73a8a0aecdd35",
    "xsections/nusigma_p_sigma_CC.dat": "b0c761caf45fbf9514774a90705a57a392f77470f741f7aece3ef3dbba332e31",
    "xsections/nusigma_p_sigma_NC.dat": "72b3ac65f03760a1386f9d48927a9227a2e031a6270331341197ce73468e52b4",
    "xsections/nusigma_n_sigma_CC.dat": "e4f71f26108e024d33df19cbe97c0fa12b18cf64f8f1d6d859f93150b5ef687c",
    "xsections/nusigma_n_sigma_NC.dat": "e1b0380f845424a97858188025f2d69dc8a52919b2688f3aa805edb2d6206314",
    "xsections/nusigma_dsde_CC.dat": "e44efc2c492df3d5f9218e68a6fdaf01d7a379d110c43a194814ee61ebf0dd3a",
    "xsections/nusigma_dsde_NC.dat": "1760d8404dcb21f27451bb6d86605ce6251574d633cf60059b6b0cbd2659c70a",

    # WCG24 nuclear cross sections (~195 MB total)
    "xsections/wcg24_base_isoscalar.h5": "c29246166a3b5f4d3200105badccf344cd8eb1629e68d71523137b208aee7b23",
    "xsections/wcg24_base_proton.h5": "90b29087c1d806afd124adf72312a91d937dccd1438db3a128b410b80af326f9",
    "xsections/wcg24_base_neutron.h5": "c982a6f5387aa001fb2dd8f290cf65e72d7127b54e15a444e2573cc68c379f44",
    "xsections/wcg24_oxygen.h5": "da041017f911eb3a426c427ee7f81422951ee8f33c6d6eabfbaece44d79fa626",
    "xsections/wcg24_carbon.h5": "9e5bf9f383acfca9aa21f2fc79967a85f7fed2d4c4f5d640f1ef67f3b1e63ee0",
    "xsections/wcg24_sodium.h5": "58806b697c320ccadfa1ae25ed26614931f72a0e01f3c48f15a654e82a6880dc",
    "xsections/wcg24_magnesium.h5": "6de920a7186feb250e689bfb509063bbc585080e5a7f684a455f93053c22c214",
    "xsections/wcg24_aluminum.h5": "29f6da7253bea53ea93fce9dafca80275f39a4c8272562d67bfc866561a89ce6",
    "xsections/wcg24_silicon.h5": "cf1d82f151cc3c38264102da4a8af97b24b0195c69bdef9bb6775be758e708a9",
    "xsections/wcg24_sulfur.h5": "cf1d82f151cc3c38264102da4a8af97b24b0195c69bdef9bb6775be758e708a9",
    "xsections/wcg24_calcium.h5": "9e5bf9f383acfca9aa21f2fc79967a85f7fed2d4c4f5d640f1ef67f3b1e63ee0",
    "xsections/wcg24_iron.h5": "da041017f911eb3a426c427ee7f81422951ee8f33c6d6eabfbaece44d79fa626",
    "xsections/wcg24_nickel.h5": "7249cfd32d8d826321f02101f449d5ffe18dd6d237bc348c4ab6a71368ecd2c7",
    "xsections/wcg24_lead.h5": "a6f31275d42f070f444f96d97fb577fc79b5986371f85aa451caf849674817dd",

    # Atmospheric production profiles (~9 MB)
    "atmos_prod/PROD_MODEL_MCEQ.dat": "9bdb755495842c155844e7651b60bef716f478e463016acd98c60581ed79bc06",
}

# Group files by category for selective downloads
FILE_GROUPS = {
    "csms": [f for f in REGISTRY if "csms" in f],
    "nusigma": [f for f in REGISTRY if "nusigma" in f],
    "wcg24": [f for f in REGISTRY if "wcg24" in f],
    "atmos": [f for f in REGISTRY if "atmos_prod" in f],
    "all": list(REGISTRY.keys()),
}


def get_data_home():
    """
    Get the path to the nuSQuIDS data cache directory.

    The cache location can be customized via the NUSQUIDS_DATA_HOME environment
    variable. By default, it uses:
    - Linux/macOS: ~/.local/share/nuSQuIDS
    - Windows: %LOCALAPPDATA%/nuSQuIDS

    Returns
    -------
    Path
        Path to the data cache directory
    """
    if "NUSQUIDS_DATA_HOME" in os.environ:
        return Path(os.environ["NUSQUIDS_DATA_HOME"])

    if sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    else:
        base = Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share"))

    return base / "nuSQuIDS"


def _create_pooch():
    """Create the pooch instance for managing data downloads."""
    if not HAS_POOCH:
        return None

    return pooch.create(
        path=get_data_home(),
        base_url=BASE_URL,
        version=DATA_VERSION,
        version_dev="main",
        registry=REGISTRY,
        env="NUSQUIDS_DATA_HOME",
    )


# Global pooch instance (lazy initialization)
_POOCH = None


def _get_pooch():
    """Get or create the global pooch instance."""
    global _POOCH
    if _POOCH is None:
        _POOCH = _create_pooch()
    return _POOCH


def fetch_file(filename, progressbar=True):
    """
    Fetch a data file, downloading it if necessary.

    Parameters
    ----------
    filename : str
        Name of the file relative to the data directory (e.g., "xsections/csms.h5")
    progressbar : bool
        Whether to show a progress bar during download

    Returns
    -------
    str
        Full path to the cached file

    Raises
    ------
    ImportError
        If pooch is not installed
    ValueError
        If the filename is not in the registry
    """
    if not HAS_POOCH:
        raise ImportError(
            "The 'pooch' package is required to download data files. "
            "Install it with: pip install pooch"
        )

    if filename not in REGISTRY:
        raise ValueError(
            f"Unknown data file: {filename}. "
            f"Available files: {list(REGISTRY.keys())}"
        )

    pup = _get_pooch()
    downloader = pooch.HTTPDownloader(progressbar=progressbar)

    return pup.fetch(filename, downloader=downloader)


def fetch_group(group="all", progressbar=True):
    """
    Fetch a group of data files.

    Parameters
    ----------
    group : str
        Name of the file group to download. Options:
        - "csms": CSMS cross sections (~50MB)
        - "nusigma": NuSigma cross sections (~30MB)
        - "wcg24": WCG24 nuclear cross sections (~200MB)
        - "atmos": Atmospheric production profiles (~5MB)
        - "all": All data files (~285MB)
    progressbar : bool
        Whether to show a progress bar during download

    Returns
    -------
    list
        List of paths to the downloaded files
    """
    if group not in FILE_GROUPS:
        raise ValueError(
            f"Unknown file group: {group}. "
            f"Available groups: {list(FILE_GROUPS.keys())}"
        )

    files = FILE_GROUPS[group]
    paths = []

    for filename in files:
        print(f"Fetching {filename}...")
        path = fetch_file(filename, progressbar=progressbar)
        paths.append(path)

    return paths


def get_xsection_path(name):
    """
    Get the path to a cross-section file, downloading if necessary.

    Parameters
    ----------
    name : str
        Name of the cross-section file (e.g., "csms.h5", "wcg24_oxygen.h5")

    Returns
    -------
    str
        Full path to the cross-section file
    """
    filename = f"xsections/{name}"
    if filename not in REGISTRY:
        # Try without the xsections/ prefix
        if name in REGISTRY:
            filename = name
        else:
            raise ValueError(f"Unknown cross-section file: {name}")

    return fetch_file(filename)


def is_data_available(filename=None):
    """
    Check if data files are available locally.

    Parameters
    ----------
    filename : str, optional
        Specific file to check. If None, checks if any cross-section files exist.

    Returns
    -------
    bool
        True if the data is available locally
    """
    data_home = get_data_home()

    if filename is not None:
        return (data_home / DATA_VERSION / filename).exists()

    # Check if any CSMS files exist (minimum required for basic operation)
    csms_files = ["xsections/csms_proton.h5", "xsections/csms_neutron.h5"]
    return all((data_home / DATA_VERSION / f).exists() for f in csms_files)


def clear_cache():
    """
    Clear all cached data files.

    This removes all downloaded data files from the cache directory.
    Use with caution as files will need to be re-downloaded.
    """
    import shutil
    data_home = get_data_home()
    if data_home.exists():
        shutil.rmtree(data_home)
        print(f"Cleared cache at {data_home}")
    else:
        print("Cache directory does not exist")


def print_cache_info():
    """Print information about the data cache."""
    data_home = get_data_home()
    print(f"Data cache location: {data_home}")
    print(f"Data version: {DATA_VERSION}")

    if not data_home.exists():
        print("Cache directory does not exist (no data downloaded)")
        return

    # Calculate total size
    total_size = 0
    file_count = 0
    for root, dirs, files in os.walk(data_home):
        for f in files:
            fp = Path(root) / f
            total_size += fp.stat().st_size
            file_count += 1

    print(f"Cached files: {file_count}")
    print(f"Total size: {total_size / 1024 / 1024:.1f} MB")


def main():
    """CLI entry point for downloading data files."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Download nuSQuIDS data files (cross sections, etc.)"
    )
    parser.add_argument(
        "group",
        nargs="?",
        default="csms",
        choices=list(FILE_GROUPS.keys()),
        help="Which data files to download (default: csms)"
    )
    parser.add_argument(
        "--info",
        action="store_true",
        help="Print cache information and exit"
    )
    parser.add_argument(
        "--clear",
        action="store_true",
        help="Clear the data cache"
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bars"
    )

    args = parser.parse_args()

    if args.info:
        print_cache_info()
        return

    if args.clear:
        clear_cache()
        return

    if not HAS_POOCH:
        print("Error: pooch is required to download data files.")
        print("Install it with: pip install pooch")
        sys.exit(1)

    print(f"Downloading {args.group} data files...")
    print(f"Cache location: {get_data_home()}")
    print()

    try:
        paths = fetch_group(args.group, progressbar=not args.no_progress)
        print()
        print(f"Successfully downloaded {len(paths)} files.")
    except Exception as e:
        print(f"Error downloading files: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
