cmake_minimum_required(VERSION 3.10)
project(TestFFTW)

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

# Add the include directories for 2decomp_C and FFTW
include_directories(
    ${FFTW_INCLUDE_DIRS}    # FFTW include path (via find_package or fallback)
    ./2decomp_C             # Include path for 2decomp_C headers

    ./
    )

# Add the source files from 2decomp_C
set(DECOMP_SRC
    ./2decomp_C/C2Decomp.cpp
    ./2decomp_C/Alloc.cpp
    ./2decomp_C/Best2DGrid.cpp
    ./2decomp_C/Halo.cpp
    ./2decomp_C/IO.cpp
    ./2decomp_C/MemSplitMerge.cpp
    ./2decomp_C/TransposeX2Y.cpp
    ./2decomp_C/TransposeY2X.cpp
    ./2decomp_C/TransposeY2Z.cpp
    ./2decomp_C/TransposeZ2Y.cpp

    ./poissonsolver.cpp
)

# Set the compiler to mpic++ for parallel compilation
set(CMAKE_CXX_COMPILER mpic++)

# Create the executable and include all .cpp files found in ./lib
add_executable(test_fft N_with_fft_3d.cpp ${DECOMP_SRC})

# Link the FFTW library and the math library
target_link_libraries(test_fft ${FFTW_LIBRARIES} m)
