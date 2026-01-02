nuSQuIDS
========

nuSQuIDS is a neutrino oscillation library that solves the neutrino evolution equations in a variety of environments (vacuum, constant density, Earth, Sun, etc.) and supports non-standard interactions.

Try it Online
-------------

No installation required! Run an interactive Jupyter notebook example:

<a href="https://mybinder.org/v2/gh/arguelles/nuSQuIDS/master?filepath=resources%2Fpython%2Fexample%2FnuSQUIDS-DEMO.ipynb" target="_blank"><img src="https://mybinder.org/badge_logo.svg" alt="Binder"></a>

Installation (Python)
---------------------

The Python bindings can be installed via pip:

```bash
pip install .
```

**Prerequisites:** GSL, HDF5, and SQuIDS must be installed on your system (see Prerequisites section below).

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
* SQuIDS (>= 1.2): https://github.com/jsalvado/SQuIDS/

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

Installing
----------

Install the library (default location: `/usr/local`):

	make install

To change the installation prefix:

	./configure --prefix=$HOME

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
