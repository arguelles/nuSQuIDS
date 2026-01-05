cmake_minimum_required(VERSION 3.15)
project(nuSQuIDS LANGUAGES CXX)

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
pkg_check_modules(SQuIDS QUIET squids)

if(SQuIDS_FOUND)
    message(STATUS "Using system SQuIDS")
    # Derive the prefix from the first include dir
    list(GET SQuIDS_INCLUDE_DIRS 0 SQuIDS_INCLUDE_DIR)
    get_filename_component(SQuIDS_PREFIX "${SQuIDS_INCLUDE_DIR}" DIRECTORY)
    set(SQuIDS_INCLUDE ${SQuIDS_INCLUDE_DIRS})
    set(SQuIDS_LIBS    ${SQuIDS_LIBRARIES})
else()
    message(STATUS "SQuIDS installation not found. Will clone and build SQuIDS")
    set(SQuIDS_PREFIX "${CMAKE_BINARY_DIR}/_deps/squids/install")

    ExternalProject_Add(squids_external
      GIT_REPOSITORY https://github.com/jsalvado/SQuIDS.git
      GIT_TAG        master
      PREFIX         ${CMAKE_BINARY_DIR}/_deps/squids
      BUILD_IN_SOURCE 1            # (switch OFF if SQuIDS supports out-of-tree)
      CONFIGURE_COMMAND
          ./configure
              --prefix=${SQuIDS_PREFIX}
              --with-gsl=${GSL_PREFIX}
      BUILD_COMMAND   make -j${CMAKE_BUILD_PARALLEL_LEVEL}
      INSTALL_COMMAND make install
      UPDATE_DISCONNECTED 1
      BUILD_ALWAYS     OFF
    )

    set(SQuIDS_INCLUDE ${SQuIDS_PREFIX}/include)
    find_library(SQuIDS_LIBS
        NAMES SQuIDS
        PATHS ${SQuIDS_PREFIX}/lib
        NO_DEFAULT_PATH)
endif()

# Extract single include path from GSL (first element of the list)
list(GET GSL_INCLUDE_DIRS 0 GSL_INCLUDE_DIR)
get_filename_component(GSL_PREFIX "${GSL_INCLUDE_DIR}" DIRECTORY)

# Extract single include path from HDF5
list(GET HDF5_INCLUDE_DIRS 0 HDF5_INCLUDE_DIR)
get_filename_component(HDF5_PREFIX "${HDF5_INCLUDE_DIR}" DIRECTORY)

# Extract single include path from SQuIDS (via pkg-config)
#list(GET SQuIDS_INCLUDE_DIRS 0 SQuIDS_INCLUDE_DIR)
#get_filename_component(SQuIDS_PREFIX "${SQuIDS_INCLUDE_DIR}" DIRECTORY)

# Directory where Autotools “make install” will land
set(NUSQUIDS_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/_deps/nusquids_core-install)

# -------------- Build core library via configure + make ---------------------

# Get python-config path that matches the Python interpreter
# This is important to avoid mismatches between python and python-config
get_filename_component(PYTHON_BIN_DIR ${Python_EXECUTABLE} DIRECTORY)
find_program(PYTHON_CONFIG
    NAMES python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}-config python${Python_VERSION_MAJOR}-config python-config
    HINTS ${PYTHON_BIN_DIR}
    NO_DEFAULT_PATH
)
if(NOT PYTHON_CONFIG)
    # Fallback to searching in PATH
    find_program(PYTHON_CONFIG
        NAMES python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}-config python${Python_VERSION_MAJOR}-config python-config
    )
endif()
message(STATUS "Using python-config: ${PYTHON_CONFIG}")

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
        --python-config=${PYTHON_CONFIG}
        --python-module-dir=${NUSQUIDS_INSTALL_PREFIX}/lib
  BUILD_COMMAND  make -j${CMAKE_BUILD_PARALLEL_LEVEL} && make python -j${CMAKE_BUILD_PARALLEL_LEVEL}
  INSTALL_COMMAND make install PREFIX=${NUSQUIDS_INSTALL_PREFIX} && make python-install PREFIX=${NUSQUIDS_INSTALL_PREFIX}
  BUILD_IN_SOURCE 1
)

if (NOT SQuIDS_FOUND)
    add_dependencies(nusquids_core squids_external)
endif()

# -------------- Installation of C++ Library ---------------------------------

# Install C++ headers into the wheel/system prefix (for downstream use)
install(DIRECTORY ${NUSQUIDS_INSTALL_PREFIX}/include/
        DESTINATION include)

# Install the compiled C++ library (.a/.so) into the wheel/system prefix
install(DIRECTORY ${NUSQUIDS_INSTALL_PREFIX}/lib/
        DESTINATION lib)

# Also install libSQuIDS to lib/ when built from source
if(NOT SQuIDS_FOUND)
    install(DIRECTORY ${SQuIDS_PREFIX}/lib/
            DESTINATION lib
            FILES_MATCHING PATTERN "libSQuIDS*")
endif()

# Install nusquids-config script
install(PROGRAMS ${NUSQUIDS_INSTALL_PREFIX}/bin/nusquids-config
        DESTINATION bin)

# -------------- Installation of Python Library ---------------------------------

# 1. Pure-Python helpers + __init__.py  →  nuSQuIDS/
install(
    DIRECTORY   ${CMAKE_SOURCE_DIR}/resources/python/bindings/
    DESTINATION nuSQuIDS
    FILES_MATCHING
        PATTERN "*.py"       # nuSQUIDSTools.py, __init__.py, etc.
)

# 2. Compiled extension built by `make python(-install)`  →  nuSQuIDS/
install(DIRECTORY ${NUSQUIDS_INSTALL_PREFIX}/lib/
        DESTINATION ${PROJECT_NAME}  # "nuSQuIDS"
        FILES_MATCHING PATTERN "*.so" PATTERN "*.dylib" PATTERN "*.pyd")

# 3. Install libSQuIDS (when built from source) into nuSQuIDS/ for bundling
#    This is needed for wheel distribution where SQuIDS isn't a system library
if(NOT SQuIDS_FOUND)
    install(DIRECTORY ${SQuIDS_PREFIX}/lib/
            DESTINATION ${PROJECT_NAME}
            FILES_MATCHING PATTERN "libSQuIDS*")
endif()

# -------------- macOS Library Path Fixups ------------------------------------
# For wheel distribution, we need to fix library paths to use @loader_path
# so the Python extension can find libnuSQuIDS and libSQuIDS in the same directory

if(APPLE)
    # Create a script that will be run after installation to fix paths
    install(CODE "
        set(NUSQUIDS_DIR \"\${CMAKE_INSTALL_PREFIX}/nuSQuIDS\")

        # Helper macro to fix library references in a binary
        macro(fix_lib_refs BINARY)
            if(EXISTS \${BINARY})
                execute_process(
                    COMMAND otool -L \${BINARY}
                    OUTPUT_VARIABLE OTOOL_OUTPUT
                )

                # Fix references to libSQuIDS (may be absolute paths from build)
                string(REGEX MATCHALL \"[^\\n]*libSQuIDS[^\\n]*\" SQUIDS_REFS \"\${OTOOL_OUTPUT}\")
                foreach(REF \${SQUIDS_REFS})
                    if(REF MATCHES \"(/[^ ]+libSQuIDS[^ ]*\\\\.dylib)\")
                        set(OLD_PATH \${CMAKE_MATCH_1})
                        get_filename_component(LIB_NAME \${OLD_PATH} NAME)
                        execute_process(
                            COMMAND install_name_tool -change
                                \${OLD_PATH}
                                @loader_path/\${LIB_NAME}
                                \${BINARY}
                            ERROR_QUIET
                        )
                    endif()
                endforeach()
                # Also fix non-absolute libSQuIDS references
                execute_process(
                    COMMAND install_name_tool -change
                        libSQuIDS.1.dylib
                        @loader_path/libSQuIDS.1.dylib
                        \${BINARY}
                    ERROR_QUIET
                )

                # Fix references to libnuSQuIDS (may be absolute paths from build)
                string(REGEX MATCHALL \"[^\\n]*libnuSQuIDS[^\\n]*\" NUSQUIDS_REFS \"\${OTOOL_OUTPUT}\")
                foreach(REF \${NUSQUIDS_REFS})
                    if(REF MATCHES \"(/[^ ]+libnuSQuIDS[^ ]*\\\\.dylib)\")
                        set(OLD_PATH \${CMAKE_MATCH_1})
                        get_filename_component(LIB_NAME \${OLD_PATH} NAME)
                        execute_process(
                            COMMAND install_name_tool -change
                                \${OLD_PATH}
                                @loader_path/\${LIB_NAME}
                                \${BINARY}
                            ERROR_QUIET
                        )
                    endif()
                endforeach()
                # Also fix non-absolute libnuSQuIDS references
                execute_process(
                    COMMAND install_name_tool -change
                        libnuSQuIDS.1.dylib
                        @loader_path/libnuSQuIDS.1.dylib
                        \${BINARY}
                    ERROR_QUIET
                )
            endif()
        endmacro()

        # Fix the Python extension nuSQuIDS.so
        set(NUSQUIDS_SO \"\${NUSQUIDS_DIR}/nuSQuIDS.so\")
        if(EXISTS \${NUSQUIDS_SO})
            message(STATUS \"Fixing library paths in nuSQuIDS.so\")
            fix_lib_refs(\${NUSQUIDS_SO})
            # Fix lib/nuSQuIDS.so reference (incorrect id path)
            execute_process(
                COMMAND install_name_tool -change
                    lib/nuSQuIDS.so
                    @loader_path/nuSQuIDS.so
                    \${NUSQUIDS_SO}
                ERROR_QUIET
            )
            # Set the install name id properly
            execute_process(
                COMMAND install_name_tool -id
                    @loader_path/nuSQuIDS.so
                    \${NUSQUIDS_SO}
                ERROR_QUIET
            )
        endif()

        # Fix all libnuSQuIDS dylibs
        file(GLOB NUSQUIDS_DYLIBS \"\${NUSQUIDS_DIR}/libnuSQuIDS*.dylib\")
        foreach(DYLIB \${NUSQUIDS_DYLIBS})
            message(STATUS \"Fixing library paths in \${DYLIB}\")
            fix_lib_refs(\${DYLIB})
            get_filename_component(DYLIB_NAME \${DYLIB} NAME)
            execute_process(
                COMMAND install_name_tool -id
                    @loader_path/\${DYLIB_NAME}
                    \${DYLIB}
                ERROR_QUIET
            )
        endforeach()

        # Fix all libSQuIDS dylibs (when bundled)
        file(GLOB SQUIDS_DYLIBS \"\${NUSQUIDS_DIR}/libSQuIDS*.dylib\")
        foreach(DYLIB \${SQUIDS_DYLIBS})
            message(STATUS \"Fixing library paths in \${DYLIB}\")
            fix_lib_refs(\${DYLIB})
            get_filename_component(DYLIB_NAME \${DYLIB} NAME)
            execute_process(
                COMMAND install_name_tool -id
                    @loader_path/\${DYLIB_NAME}
                    \${DYLIB}
                ERROR_QUIET
            )
        endforeach()
    ")
endif()

# -------------- Linux Library Path Fixups ------------------------------------
# For wheel distribution on Linux, we need to set RPATH to $ORIGIN
# so the Python extension can find shared libraries in the same directory

if(UNIX AND NOT APPLE)
    install(CODE "
        set(NUSQUIDS_DIR \"\${CMAKE_INSTALL_PREFIX}/nuSQuIDS\")

        # Fix RPATH for nuSQuIDS.so
        set(NUSQUIDS_SO \"\${NUSQUIDS_DIR}/nuSQuIDS.so\")
        if(EXISTS \${NUSQUIDS_SO})
            message(STATUS \"Setting RPATH in nuSQuIDS.so\")
            execute_process(
                COMMAND patchelf --set-rpath \\$ORIGIN \${NUSQUIDS_SO}
                ERROR_QUIET
            )
        endif()

        # Fix RPATH for all .so files in the nuSQuIDS directory
        file(GLOB ALL_SO_FILES \"\${NUSQUIDS_DIR}/*.so\" \"\${NUSQUIDS_DIR}/*.so.*\")
        foreach(SO_FILE \${ALL_SO_FILES})
            message(STATUS \"Setting RPATH in \${SO_FILE}\")
            execute_process(
                COMMAND patchelf --set-rpath \\$ORIGIN \${SO_FILE}
                ERROR_QUIET
            )
        endforeach()
    ")
endif()

