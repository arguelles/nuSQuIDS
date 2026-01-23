cmake_minimum_required(VERSION 3.15)
project(nuSQuIDSjl CXX)

# Use C++17 for compatibility with libcxxwrap-julia
# Add compatibility for std::result_of which was removed in C++20
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add compatibility macro for std::result_of (deprecated in C++17, removed in C++20)
# std::result_of<F(Args...)>::type -> std::invoke_result_t<F, Args...>
add_compile_definitions(_LIBCPP_ENABLE_CXX20_REMOVED_TYPE_TRAITS)

# Find required packages
find_package(JlCxx REQUIRED)
find_package(PkgConfig REQUIRED)

# Find nuSQuIDS dependencies
pkg_check_modules(GSL REQUIRED gsl)
pkg_check_modules(SQUIDS REQUIRED squids)

# Use local nuSQuIDS source headers (for latest C++17 compatibility fixes)
# and system library for linking
set(NUSQUIDS_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../../include")
set(NUSQUIDS_LIBRARY_DIRS "/usr/local/lib")
set(NUSQUIDS_LIBRARIES "nuSQuIDS")

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX HL)

# Create the Julia module library
add_library(nuSQuIDSjl SHARED
    src/nuSQuIDSjl.cpp
)

target_include_directories(nuSQuIDSjl PRIVATE
    ${JlCxx_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
    ${SQUIDS_INCLUDE_DIRS}
    ${NUSQUIDS_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

target_link_directories(nuSQuIDSjl PRIVATE
    ${GSL_LIBRARY_DIRS}
    ${SQUIDS_LIBRARY_DIRS}
    ${NUSQUIDS_LIBRARY_DIRS}
    ${HDF5_LIBRARY_DIRS}
)

target_link_libraries(nuSQuIDSjl
    JlCxx::cxxwrap_julia
    JlCxx::cxxwrap_julia_stl
    ${NUSQUIDS_LIBRARIES}
    ${SQUIDS_LIBRARIES}
    ${GSL_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
)

# Set output directory
set_target_properties(nuSQuIDSjl PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
    PREFIX ""
)

# Platform-specific settings
if(APPLE)
    set_target_properties(nuSQuIDSjl PROPERTIES
        SUFFIX ".dylib"
        INSTALL_RPATH "@loader_path"
    )
else()
    set_target_properties(nuSQuIDSjl PROPERTIES
        SUFFIX ".so"
        INSTALL_RPATH "$ORIGIN"
    )
endif()

# Installation
install(TARGETS nuSQuIDSjl
    LIBRARY DESTINATION lib
)
