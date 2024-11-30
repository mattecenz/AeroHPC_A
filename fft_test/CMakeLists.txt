cmake_minimum_required(VERSION 3.10)
project(TestFFTW)

# Specify the include directories
include_directories(
    /usr/local/include      # FFTW includes
    ./2decomp_C             # Include path for 2decomp_C headers
)

link_directories(/usr/local/lib)

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
)

# Set the compiler to mpic++ for parallel compilation
set(CMAKE_CXX_COMPILER mpic++)

# Create the executable and include all .cpp files found in ./lib
# add_executable(test_fft with_fft.cpp ${DECOMP_SRC})

add_executable(test_fft N_test_fftw_vector.cpp ${DECOMP_SRC})

# Link the FFTW library and the math library
target_link_libraries(test_fft fftw3 m)
