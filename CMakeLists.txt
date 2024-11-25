cmake_minimum_required(VERSION 3.20)
project(AeroHPC_A)

# C++ Standard
set(CMAKE_CXX_STANDARD 20)

# Include header files
include_directories(include)

# Include source files
file(GLOB SOURCES "src/*.cpp")

# Main file name
set(MAIN_FILE "main.cpp")

# Target file name
set(TARGET_FILE "AeroHPC_A")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -O3")

add_executable("${TARGET_FILE}" "${MAIN_FILE}" ${SOURCES})

#ForcingTerm FLAG
add_definitions(-DForcingT)


#Logging FLAG
add_definitions(-DLogging)
