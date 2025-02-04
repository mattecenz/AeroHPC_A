cmake_minimum_required(VERSION 3.20)
project(Div_test)

# C++ Standard
set(CMAKE_CXX_STANDARD 23)

# Include header files
include_directories(../../include)

set(MAIN_FILE "main.cpp")

# Target file name
set(TARGET_FILE "Div_test")

add_executable("${TARGET_FILE}" "${MAIN_FILE}" "../../src/Traits.cpp" "../../src/GridData.cpp")