i3_project(nuSQuIDS CXX)

## build our library
i3_add_library(nuSQuIDS
  private/nuSQuIDS/*.cpp
  USE_PROJECTS SQuIDS
  USE_TOOLS gsl hdf5
  WITHOUT_I3_HEADERS)
set_source_files_properties(${PROJECT_SOURCE_DIR}/private/${PROJECT_NAME}/body.cpp     COMPILE_FLAGS "-std=c++11")
set_source_files_properties(${PROJECT_SOURCE_DIR}/private/${PROJECT_NAME}/nuSQUIDS.cpp COMPILE_FLAGS "-std=c++11")

## build the pybindings
i3_add_pybindings(nuSQuIDS
  private/pybindings/nuSQUIDSpy.cpp
  USE_PROJECTS SQuIDS nuSQuIDS
  USE_TOOLS gsl hdf5 python boostnumpy)
## i3_add_pybindings() doesn't support eh WITHOUT_I3_HEADERS flag, so this disables it by hand
set_target_properties(${PROJECT_NAME}-pybindings
  PROPERTIES
  COMPILE_FLAGS "")
set_source_files_properties(${PROJECT_SOURCE_DIR}/private/pybindings/nuSQUIDSpy.cpp COMPILE_FLAGS "-std=c++11")

## define the tests
enable_testing()
file(GLOB tests ${PROJECT_SOURCE_DIR}/private/test/*.test.cpp)
foreach(t ${tests})
  get_filename_component(f ${t} NAME_WE)
  i3_executable("${f}"
    "${t}"
    USE_PROJECTS SQuIDS nuSQuIDS
    USE_TOOLS gsl hdf5
    WITHOUT_I3_HEADERS)
  set_source_files_properties(${t} COMPILE_FLAGS "-std=c++11")
  add_dependencies(test-bins ${PROJECT_NAME}-${f})
  add_test(NAME "${PROJECT_NAME}::${f}"
           WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
	   COMMAND "${PROJECT_NAME}-${f}")
endforeach()
