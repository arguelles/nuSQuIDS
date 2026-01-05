nuSQuIDS
========

nuSQuIDS is a neutrino oscillation library that solves the neutrino evolution equations in a variety of environments (vacuum, constant density, Earth, Sun, etc.) and supports non-standard interactions.

Try it Online
-------------

No installation required! Run an interactive Jupyter notebook example:

<a href="https://mybinder.org/v2/gh/arguelles/nuSQuIDS/master?filepath=resources%2Fpython%2Fexample%2FnuSQUIDS-DEMO.ipynb" target="_blank"><img src="https://mybinder.org/badge_logo.svg" alt="Binder"></a>

Installation (Python)
---------------------

### From PyPI (Recommended)

The easiest way to install nuSQuIDS is from PyPI:

```bash
pip install nusquids
```

This installs pre-built wheels with all dependencies bundled - no need to install GSL, HDF5, or SQuIDS separately.

**Supported platforms:**
- Linux (x86_64, aarch64) - Python 3.9+
- macOS (arm64) - Python 3.10+

### From Source

To build from source (e.g., for development or unsupported platforms):

```bash
pip install .
```

**Prerequisites for source builds:** GSL, HDF5, and SQuIDS must be installed on your system (see Prerequisites section below).

After installation:
```python
import nuSQuIDS as nsq
print("nuSQuIDS imported successfully!")
```

Prerequisites
-------------

The following system libraries are required:

* GSL (>= 1.15): http://www.gnu.org/software/gsl/
* HDF5 with C bindings: http://www.hdfgroup.org/HDF5/
* SQuIDS (> 1.2): https://github.com/jsalvado/SQuIDS/

For Python bindings, you also need:
* numpy: http://www.numpy.org/

Optional (recommended for interactive use):
* matplotlib: http://matplotlib.org/
* ipython: http://ipython.org/

Documentation
-------------

The LaTeX and PDF user manual is included in the `doc/` folder.
To generate documentation from source: `make docs`

---

Installation from Source (C++ Library)
======================================

The following instructions are for building the C++ library from source.
**If you only need the Python interface, `pip install .` is sufficient.**

Configuration
-------------

The path for the GSL libraries can be specified by running:

	./configure --with-gsl-incdir=GSL_include_path --with-gsl-libdir=GSL_library_path

or, assuming a standard installation into 'include' and 'lib' subdirectories within a common prefix:

	./configure --with-gsl=GSL_prefix

The path for the HDF5 libraries can be specified by running:

	./configure --with-hdf5-incdir=HDF5_include_path --with-hdf5-libdir=HDF5_library_path

or more simply:

	./configure --with-hdf5=HDF5_prefix

Finally, the path for SQuIDS can be specified by:

	./configure --with-squids-incdir=SQuIDS_include_path --with-squids-libdir=SQuIDS_library_path

or commonly just:

	./configure --with-squids=SQuIDS_prefix

Building
--------

Once configuration is complete, compile the library:

	make

Run unit tests to verify correct behavior:

	make test

Compile example programs:

	make examples

The resulting example executables can be found in the various subdirectories of `examples/`.

Examples
--------

The `examples/` directory contains both C++ and Python versions of example programs:

| Directory | Description |
|-----------|-------------|
| `Single_energy/` | Single energy mode oscillations |
| `Multiple_energy/` | Multiple energy mode with power-law spectrum |
| `Bodies/` | Various body types (Earth, EarthAtm, Vacuum, etc.) |
| `Constant_density_layers/` | Multi-layer propagation |
| `Atm_default/` | Atmospheric neutrino propagation (nuSQUIDSAtm) |
| `HDF5_Write_Read/` | State serialization to HDF5 |
| `Astrophysical_neutrino_flavor_ratio/` | Astrophysical flavor ratio calculation |
| `NSI/` | Non-standard interactions (C++ only) |
| `LV/` | Lorentz violation (C++ only) |
| `Xsections/` | Custom cross-sections (C++ only) |

Run C++ examples after building:

	./examples/Single_energy/single_energy

Run Python examples:

	python examples/Single_energy/main.py

Benchmarks
----------

nuSQuIDS includes a benchmark suite to measure performance:

	make benchmark       # Full benchmark (more iterations)
	make benchmark-quick # Quick benchmark

The benchmark tests single energy, multiple energy, and atmospheric modes with various physics configurations (with/without interactions, Glashow resonance, tau regeneration).

Installing
----------

Install the library (default location: `/usr/local`):

	make install

To change the installation prefix:

	./configure --prefix=$HOME

This installs:
- Headers to `$PREFIX/include/nuSQuIDS/`
- Libraries to `$PREFIX/lib/`
- The `nusquids-config` script to `$PREFIX/bin/`
- Physics data to `$PREFIX/share/nuSQuIDS/`

Using nusquids-config
---------------------

After installation, use `nusquids-config` to get compiler and linker flags for your own programs:

	# Get compiler flags
	nusquids-config --cflags

	# Get linker flags
	nusquids-config --libs

	# Compile your own program
	g++ $(nusquids-config --cflags) -c myprogram.cpp
	g++ myprogram.o $(nusquids-config --libs) -o myprogram

	# Or in one step
	g++ $(nusquids-config --cflags) myprogram.cpp $(nusquids-config --libs) -o myprogram

Available options: `--prefix`, `--incdir`, `--libdir`, `--cflags`, `--libs`, `--version`

Python Bindings (Manual Build)
------------------------------

To build Python bindings manually (instead of using `pip install .`), pass the `--with-python-bindings` option to configure:

	./configure --with-python-bindings

By default this uses pybind11. To use Boost.Python instead:

	./configure --with-boost-python-bindings

Additional options:

	# Specify pybind11 headers location
	./configure --with-pybind-incdir=location_of_pybind11

	# Specify Boost location (for Boost.Python)
	./configure --with-boost=boost_prefix

	# Specify Python executable
	./configure --python-bin=PYTHON_EXECUTABLE

	# Specify python-config for compilation flags
	./configure --python-config=PYTHON_CONFIG_EXECUTABLE

	# Specify installation directory for Python module
	./configure --python-module-dir=PYTHON_MODULE_INSTALLATION

Build and install the Python bindings:

	make python
	make python-install

Alternatively, add the bindings to your PYTHONPATH:

	export PYTHONPATH=$(PATH_TO_nuSQUIDS)/resources/python/bindings/:$PYTHONPATH

---

Citation
--------

If you use nuSQuIDS in your research, please cite our paper:

```bibtex
@article{Arguelles:2021twb,
    author = {Arg{\"u}elles, Carlos A. and Salvado, Jordi and Weaver, Christopher N.},
    title = "{nuSQuIDS: A toolbox for neutrino propagation}",
    eprint = "2112.13804",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.cpc.2022.108346",
    journal = "Comput. Phys. Commun.",
    volume = "277",
    pages = "108346",
    year = "2022"
}
```

nuSQuIDS is built on top of the SQuIDS library. If appropriate, please also consider citing the SQuIDS papers:

```bibtex
@article{Arguelles:2020hss,
    author = {Arg{\"u}elles, Carlos A. and Salvado, Jordi and Weaver, Christopher N.},
    title = "{A Simple Quantum Integro-Differential Solver (SQuIDS)}",
    doi = "10.1016/j.cpc.2020.107405",
    journal = "Comput. Phys. Commun.",
    volume = "255",
    pages = "107405",
    year = "2020"
}

@article{Delgado:2014lyt,
    author = "Arguelles, Carlos and Salvado, Jordi and Weaver, Christopher N.",
    title = "{A Simple Quantum Integro-Differential Solver (SQuIDS)}",
    eprint = "1412.3832",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.cpc.2015.06.022",
    journal = "Comput. Phys. Commun.",
    volume = "196",
    pages = "569--591",
    year = "2015"
}
```
