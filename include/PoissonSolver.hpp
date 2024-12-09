#ifndef AEROHPC_A_FAST_POISSON_SOLVER
#define AEROHPC_A_FAST_POISSON_SOLVER

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <fftw3.h>
#include <chrono>

#include "C2Decomp.hpp"

class poissonSolver {
private:
    int N;                  // Grid size
    double L;               // Domain length
    double dx;              // Grid spacing
    double *b;              // Input vector (right-hand side of Ax=b)
    double *u1, *u2, *u3;   // Intermediate buffers
    C2Decomp *c2d;          // Decomposition object
    double xSize[3], ySize[3], zSize[3];
    bool periodicBC[3];
    int pRow, pCol;
    fftw_plan fft_plan_x, fft_plan_y, fft_plan_z;
    fftw_plan ifft_plan_x, ifft_plan_y, ifft_plan_z;

    void initializeGrid();
    void performFFT();
    void performIFFT();
    void solveEigenvalues();

public:
    poissonSolver(int N, double L, double *b, bool periodicBC[3], int pRow = 0, int pCol = 0);
    ~poissonSolver();

    void solve(Real *X); // Updates the input array `X` with the solution
};

#endif
