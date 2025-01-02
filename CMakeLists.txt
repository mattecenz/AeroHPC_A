cmake_minimum_required(VERSION 3.20)
project(AeroHPC_A)

# C++ Standard
set(CMAKE_CXX_STANDARD 23)

# Optionally set default include and library paths (fallback)
set(DEFAULT_FFTW_INCLUDE_PATH /usr/local/include)
set(DEFAULT_FFTW_LIBRARY_PATH /usr/local/lib)

# Try to find FFTW using find_package
find_package(FFTW QUIET COMPONENTS SINGLE DOUBLE)

# Check if FFTW was found, otherwise use the fallback paths
if(FFTW_FOUND)
    message(STATUS "FFTW found via find_package.")
    set(FFTW_LIBRARIES FFTW::FFTW3 FFTW::FFTW3F)
    set(FFTW_INCLUDE_DIRS "") # Already handled by FFTW imported targets
else()
    message(WARNING "FFTW not found via find_package. Falling back to default paths.")
    set(FFTW_INCLUDE_DIRS ${DEFAULT_FFTW_INCLUDE_PATH})
    set(FFTW_LIBRARIES ${DEFAULT_FFTW_LIBRARY_PATH}/libfftw3.a)
endif()

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

# Add the include directories for FFTW
include_directories(${FFTW_INCLUDE_DIRS})

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

# Link the FFTW library
target_link_libraries(${TARGET_FILE} ${FFTW_LIBRARIES})