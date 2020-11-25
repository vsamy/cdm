# Version minimum
cmake_minimum_required(VERSION 3.8.2) # For c++17

set(CMAKE_CXX_STANDARD 17)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(PROJECT_NAME cdm)
set(PROJECT_DESCRIPTION "JRL cdm (Comprehensive Dynamics Model) based on coma library")
set(PROJECT_URL "...")
set(PROJECT_VERSION 0.1.0)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(DOXYGEN_USE_MATHJAX "YES")

include(cmake/base.cmake)
include(cmake/msvc-specific.cmake)

SET(CXX_DISABLE_WERROR True)
project(${PROJECT_NAME} CXX)

option(BUILD_TESTING "Build unit tests." ON)

if(MSVC)
    set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif()

# External dependencies
find_package(Eigen3 REQUIRED)
find_package(coma 1.0.0 REQUIRED)

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
add_subdirectory(src)

if (${BUILD_TESTING})
    add_subdirectory(test)
endif()
