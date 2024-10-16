# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Iinclude -Wall -Wextra -std=c++17

# Directories
INCLUDE_DIR = include
SRC_DIR = src
OBJ_DIR = obj

# Target executable
TARGET = main.out

# Source and object files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Default target
all: $(TARGET)

# Link object files to create the final executable
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up generated files
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

# Phony targets
.PHONY: all clean
