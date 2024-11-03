cmake_minimum_required(VERSION 3.20)
project(AeroHPC_A)

# C++ Standard
set(CMAKE_CXX_STANDARD 23)

# Include header files
include_directories(include)

# Include source files
file(GLOB SOURCES "src/*.cpp")

# Main file name
set(MAIN_FILE "main.cpp")

# Target file name
set(TARGET_FILE "AeroHPC_A")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -funroll-all-loops -fopt-info-all -march=native -mtune=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -mtune=native -fopenmp-simd -fopt-info-all\
         -funroll-loops  -flto  -flto=auto -fno-signed-zeros -fno-trapping-math -ffinite-math-only -ffast-math -fpermissive \
         -faggressive-loop-optimizations -funsafe-math-optimizations -ftree-vectorizer-verbose=2 ")

add_executable("${TARGET_FILE}" "${MAIN_FILE}" ${SOURCES})

#ForcingTerm FLAG
add_definitions(-DForcingT)
