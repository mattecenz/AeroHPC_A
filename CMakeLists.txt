cmake_minimum_required(VERSION 3.20)
project(AeroHPC_A)

# C++ Standard
set(CMAKE_CXX_STANDARD 23)

# Find the MPI compiler
find_program(MPICXX_EXECUTABLE mpic++)

if (MPICXX_EXECUTABLE)
    message(STATUS "Found mpic++ at: ${MPICXX_EXECUTABLE}")
    # Set mpic++ as the default C++ compiler
    set(CMAKE_CXX_COMPILER ${MPICXX_EXECUTABLE} CACHE FILEPATH "Path to the mpic++ compiler" FORCE)
else ()
    message(FATAL_ERROR "mpic++ not found. Please ensure MPI is installed and mpic++ is in your PATH.")
endif ()

# Example usage message
message(STATUS "Using C++ compiler: ${CMAKE_CXX_COMPILER}")

# Include header files
include_directories(include)

# Include source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Main file name
set(MAIN_FILE "main.cpp")

# Target file name
set(TARGET_FILE "AeroHPC_A")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto -Ofast -funroll-all-loops -fopenmp-simd -ffast-math ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -mtune=native -fopenmp-simd -fopt-info-all\
#         -funroll-loops  -flto  -flto=auto -fno-signed-zeros -fno-trapping-math -ffinite-math-only -ffast-math -fpermissive \
#         -faggressive-loop-optimizations -funsafe-math-optimizations -ftree-vectorizer-verbose=2")

add_executable("${TARGET_FILE}" "${MAIN_FILE}" ${SOURCES})

#ForcingTerm FLAG
add_definitions(-DForcingT)
