cmake_minimum_required(VERSION 3.15)
project(nusquids LANGUAGES CXX)

# -------------- Dependency discovery ----------------------------------------
include(ExternalProject)
find_package(PkgConfig        REQUIRED)
find_package(GSL              REQUIRED)                   # GSL::gsl  + GSL::gslcblas
find_package(HDF5 REQUIRED COMPONENTS CXX HL)             # HDF5::CXX + HDF5::HL
find_package(Python REQUIRED COMPONENTS Interpreter)      # gives PYTHON_EXECUTABLE

# Get pybind11 CMake path from Python
execute_process(
    COMMAND ${Python_EXECUTABLE} -m pybind11 --cmakedir
    OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Append to CMAKE_PREFIX_PATH so find_package can use it
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_CMAKE_DIR}")

find_package(pybind11 CONFIG REQUIRED)                    # gives pybind11::module

# SQuIDS via pkg-config (make sure squids.pc is on PKG_CONFIG_PATH)
pkg_check_modules(SQuIDS REQUIRED squids)

# Extract single include path from GSL (first element of the list)
list(GET GSL_INCLUDE_DIRS 0 GSL_INCLUDE_DIR)
get_filename_component(GSL_PREFIX "${GSL_INCLUDE_DIR}" DIRECTORY)

# Extract single include path from HDF5
list(GET HDF5_INCLUDE_DIRS 0 HDF5_INCLUDE_DIR)
get_filename_component(HDF5_PREFIX "${HDF5_INCLUDE_DIR}" DIRECTORY)

# Extract single include path from SQuIDS (via pkg-config)
list(GET SQuIDS_INCLUDE_DIRS 0 SQuIDS_INCLUDE_DIR)
get_filename_component(SQuIDS_PREFIX "${SQuIDS_INCLUDE_DIR}" DIRECTORY)

# Directory where Autotools “make install” will land
set(NUSQUIDS_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/_deps/nusquids_core-install)

# -------------- Build core library via configure + make ---------------------
ExternalProject_Add(nusquids_core
  SOURCE_DIR     ${CMAKE_SOURCE_DIR}
  CONFIGURE_COMMAND
      ${CMAKE_SOURCE_DIR}/configure
        --prefix=${NUSQUIDS_INSTALL_PREFIX}
        --with-hdf5=${HDF5_PREFIX}
        --with-gsl=${GSL_PREFIX}
        --with-squids=${SQuIDS_PREFIX}
        --with-python-bindings
        --python-bin=${Python_EXECUTABLE}
	--python-module-dir=${NUSQUIDS_INSTALL_PREFIX}/lib
  BUILD_COMMAND  make -j${CMAKE_BUILD_PARALLEL_LEVEL} && make python -j${CMAKE_BUILD_PARALLEL_LEVEL}
  INSTALL_COMMAND make install PREFIX=${NUSQUIDS_INSTALL_PREFIX} && make python-install PREFIX=${NUSQUIDS_INSTALL_PREFIX}
  BUILD_IN_SOURCE 1
)

# -------------- Installation of C++ Library ---------------------------------

# Install C++ headers into the wheel/system prefix (for downstream use)
install(DIRECTORY ${NUSQUIDS_INSTALL_PREFIX}/include/
        DESTINATION include)

# Install the compiled C++ library (.a/.so) into the wheel/system prefix
install(DIRECTORY ${NUSQUIDS_INSTALL_PREFIX}/lib/
        DESTINATION lib)

# -------------- Installation of Python Library ---------------------------------

# 1. Pure-Python helpers + __init__.py  →  nusquids/
install(
    DIRECTORY   ${CMAKE_SOURCE_DIR}/resources/python/bindings/
    DESTINATION nusquids
    FILES_MATCHING
        PATTERN "*.py"       # nuSQUIDSTools.py, __init__.py, etc.
)

# 2. Compiled extension built by `make python(-install)`  →  nusquids/
install(DIRECTORY ${NUSQUIDS_INSTALL_PREFIX}/lib/
        DESTINATION ${PROJECT_NAME}  # usually "nusquids"
        FILES_MATCHING PATTERN "*.so" PATTERN "*.dylib" PATTERN "*.pyd")
