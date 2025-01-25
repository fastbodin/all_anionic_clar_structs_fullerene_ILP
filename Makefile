# Compiler
CXX = g++-14

# Compiler flags
GUROBI = -I/Library/gurobi1103/macos_universal2/include -L/Library/gurobi1103/macos_universal2/lib -lgurobi_c++ -lgurobi110 -lm
INCLUDE = -I include/
CXX_FLAGS = -std=c++14 -O3

# Source file
SRCS = src/*.cpp

# Executable name
EXEC = build/comp_all_anionic_clar_structs

# Build rule
$(EXEC): $(SRCS)
	$(CXX) $(GUROBI) $(INCLUDE) $(CXX_FLAGS) $(SRCS) -o $(EXEC)

# Clean rule
clean:
	rm -f $(EXEC)
