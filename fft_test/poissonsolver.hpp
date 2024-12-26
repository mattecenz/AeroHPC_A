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
    double *fft_inputX, *fft_inputY, *fft_inputZ;         // FFT input buffer
    double *fft_outputX, *fft_outputY, *fft_outputZ;      // FFT input buffer

    void initializeGrid();
    void performFFT();
    void performIFFT();
    void solveEigenvalues();

public:
    poissonSolver(int N, double L, C2Decomp *c2d);
    ~poissonSolver();

    void setB(double *b);   // Setter function to update `b`, the rhs
    void solve(double *X); 
    // later on change this to Real *X
    // Updates the input array `X` with the solution
};

#endif