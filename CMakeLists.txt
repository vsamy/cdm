# Version minimum
cmake_minimum_required(VERSION 3.8.2) # For c++17

set(CMAKE_CXX_STANDARD 17)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(PROJECT_NAME cdm)
set(PROJECT_DESCRIPTION "JRL cdm (Comprehensive Dynamics Model) based on coma library")
set(PROJECT_URL "...")
set(PROJECT_VERSION 0.2.0)
# set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(DOXYGEN_USE_MATHJAX "YES")

include(cmake/base.cmake)
include(cmake/msvc-specific.cmake)

# SET(CXX_DISABLE_WERROR True)
project(${PROJECT_NAME} CXX)

option(BUILD_TESTING "Build unit tests." ON)
option(BUILD_V0 "Build the v0 library." OFF)

if(MSVC)
    set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif()

# External dependencies
add_project_dependency(coma REQUIRED)

# set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
# set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

if(NOT ${BUILD_V0})
    if(${PROJECT_VERSION} VERSION_LESS "1.0.0")
        message(WARNING "Current version is still a WIP state and is not yet usable. You can build working algorithm in the algo_v0 folder and build it by passing -DBUILD_V0=ON to cmake.")
    endif()
    add_subdirectory(src)

    if (${BUILD_TESTING})
        add_subdirectory(test)
    endif()
else()
    add_subdirectory(algo_v0)
endif()