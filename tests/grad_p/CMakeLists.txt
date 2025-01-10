cmake_minimum_required(VERSION 3.20)
project(Grad_p_test)

# C++ Standard
set(CMAKE_CXX_STANDARD 23)

# Include header fpiles
include_directories(../../include)

set(MAIN_FILE "main.cpp")

# Target file name
set(TARGET_FILE "Grad_p_test")

add_executable("${TARGET_FILE}" "${MAIN_FILE}" "../../src/Traits.cpp" "../../src/GridData.cpp")