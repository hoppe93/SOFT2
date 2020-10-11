cmake_minimum_required(VERSION 3.9..3.14)
# Required features/version
# >= 3.9
#    + Simple linking/compilation with OpenMP

project(soft VERSION 1.0
             DESCRIPTION "The Synchrotron-detecting Orbit Following Toolkit"
             LANGUAGES CXX)

option(BUILD_TESTS "Build unit tests" OFF)
option(COLOR_TERMINAL "Allow colorful SOFT output" ON)
option(OPTIMIZE_NATIVE "Enable native optimizations (non-debug mode only)" ON)
option(PROFILING "Enable gprof profiling hooks" OFF)
option(WITH_MPI "Enable MPI support" OFF)

# Add CMake modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

configure_file(
	"${PROJECT_SOURCE_DIR}/include/config.h.in"
	"${PROJECT_BINARY_DIR}/include/config.h"
)

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if (EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory")
endif()

add_subdirectory(src)
if (BUILD_TESTS)
    add_subdirectory(tests)
endif()

