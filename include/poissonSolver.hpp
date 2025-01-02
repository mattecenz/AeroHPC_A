#ifndef AEROHPC_A_FAST_POISSON_SOLVER
#define AEROHPC_A_FAST_POISSON_SOLVER

#include <fftw3.h>
#include "C2Decomp.hpp"
#include "Traits.hpp"

class poissonSolver {

// Rename the fftw library function based on Real type:
// "fftw_" functions are for double type
// "fftwf_" functions are for float type
// The renamed format is "fftwr_" (where r stands for Real)
#if REAL_USE_FLOAT
#define fftwr_malloc(a) real_p(fftwf_malloc(a))
#define fftwr_free fftwf_free
#define fftwr_plan_r2r_1d fftwf_plan_r2r_1d
#define fftwr_execute fftwf_execute
#else
#define fftwr_malloc(a) real_p(fftw_malloc(a))
#define fftwr_free fftw_free
#define fftwr_plan_r2r_1d fftw_plan_r2r_1d
#define fftwr_execute fftw_execute
#endif

private:
    int N; // Grid size
    Real L; // Domain length
    Real dx; // Grid spacing
    Real *b; // Input vector (right-hand side of Ax=b)
    Real *u1, *u2, *u3; // Intermediate buffers
    C2Decomp *c2d; // Decomposition object
    Real xSize[3], ySize[3], zSize[3];
    bool periodicBC[3];
    int pRow, pCol;
    Real *fft_inputX, *fft_inputY, *fft_inputZ; // FFT input buffer
    Real *fft_outputX, *fft_outputY, *fft_outputZ; // FFT input buffer

    void performFFT();

    void performIFFT();

    void solveEigenvalues();

public:
    poissonSolver(int N, Real L, C2Decomp *c2d);

    ~poissonSolver();

    void setB(Real *b); // Setter function to update `b`, the rhs
    void solve(Real *X);

    // later on change this to Real *X
    // Updates the input array `X` with the solution
};

#endif
