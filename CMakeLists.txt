cmake_minimum_required(VERSION 3.20)
project(AeroHPC_A)

# C++ Standard
set(CMAKE_CXX_STANDARD 23)

# Optionally set default include and library paths (fallback)
if (NOT DEFINED FFTW_INCLUDE_PATH)
    message(STATUS "Default fftw include path: /usr/local/include")
    set(DEFAULT_FFTW_INCLUDE_PATH /usr/local/include)
else()
    message(STATUS "Custom fftw include path: ${FFTW_INCLUDE_PATH}")
    set(DEFAULT_FFTW_INCLUDE_PATH ${FFTW_INCLUDE_PATH})
endif ()

if (NOT DEFINED FFTW_LIBRARY_PATH)
    message(STATUS "Default fftw library path: /usr/local/lib")
    set(DEFAULT_FFTW_LIBRARY_PATH /usr/local/lib)
else()
    message(STATUS "Custom fftw library path: ${FFTW_LIBRARY_PATH}")
    set(DEFAULT_FFTW_LIBRARY_PATH ${FFTW_LIBRARY_PATH})
endif()

# Enable use of Float values instead of Doubles
if (DEFINED USE_FLOAT)
    add_compile_definitions(REAL_USE_FLOAT=${USE_FLOAT})
endif ()

# Try to find FFTW using find_package
find_package(FFTW QUIET COMPONENTS SINGLE DOUBLE)

# Check if FFTW was found, otherwise use the fallback paths
if (FFTW_FOUND)
    message(STATUS "FFTW found via find_package.")
    set(FFTW_INCLUDE_DIRS "") # Already handled by FFTW imported targets
    if(USE_FLOAT)
        set(FFTW_LIBRARIES FFTW::FFTW3F)
    else ()
        set(FFTW_LIBRARIES FFTW::FFTW3)
    endif ()
else ()
    message(WARNING "FFTW not found via find_package. Falling back to default paths.")
    set(FFTW_INCLUDE_DIRS ${DEFAULT_FFTW_INCLUDE_PATH})
    if(USE_FLOAT)
        set(FFTW_LIBRARIES ${DEFAULT_FFTW_LIBRARY_PATH}/libfftw3f.so)
    else ()
        set(FFTW_LIBRARIES ${DEFAULT_FFTW_LIBRARY_PATH}/libfftw3.so)
    endif ()
endif ()
message(STATUS "Using fftw .so : ${FFTW_LIBRARIES}")

# Find the MPI compiler
find_program(MPICXX_EXECUTABLE mpic++)

if (MPICXX_EXECUTABLE)
    message(STATUS "Found mpic++ at: ${MPICXX_EXECUTABLE}")
    # Set mpic++ as the default C++ compiler
    set(CMAKE_CXX_COMPILER ${MPICXX_EXECUTABLE} CACHE FILEPATH "Path to the mpic++ compiler" FORCE)
else ()
    message(FATAL_ERROR "mpic++ not found. Please ensure MPI is installed and mpic++ is in your PATH.")
endif ()

message(STATUS "Using C++ compiler: ${CMAKE_CXX_COMPILER}")

# Include header files
include_directories(include)

# Add the include directories for FFTW
include_directories(${FFTW_INCLUDE_DIRS})

# Include source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Main file name
if (DEFINED MAKE_TEST)
    add_compile_definitions(TEST=1)
    add_compile_definitions(ForcingT=1)

    if(DEFINED DISABLE_PRESSURE)
        add_compile_definitions(DISABLE_PRESSURE=${DISABLE_PRESSURE})
        message(STATUS "Pressure disabled: ${DISABLE_PRESSURE}")
    endif ()
    if(DEFINED DEBUG_PRINT_BUFFERS)
        add_compile_definitions(DEBUG_PRINT_BUFFERS=${DEBUG_PRINT_BUFFERS})
        message(STATUS "Debug bugger enabled: ${DEBUG_PRINT_BUFFERS}")
    endif ()
    message(STATUS "main file: test.cpp")
    set(MAIN_FILE "test.cpp")
else ()
    message(STATUS "main file: main.cpp")
    set(MAIN_FILE "main.cpp")
endif ()

# Target file name
set(TARGET_FILE "AeroHPC_A")

if (DEFINED AGGRESSIVE_OPTIMIZATION)
    message(STATUS "Aggressive optimizations activated")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -mtune=native -fopenmp-simd -fopt-info-all\
         -funroll-loops  -flto  -flto=auto -fno-signed-zeros -fno-trapping-math -ffinite-math-only -ffast-math -fpermissive \
         -faggressive-loop-optimizations -funsafe-math-optimizations -ftree-vectorizer-verbose=2")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto -Ofast -funroll-all-loops -fopenmp-simd -ffast-math ")
endif ()

add_executable("${TARGET_FILE}" "${MAIN_FILE}" ${SOURCES})

# Link the FFTW library
target_link_libraries(${TARGET_FILE} ${FFTW_LIBRARIES})