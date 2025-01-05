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
    index_t Nx, Ny, Nz; // Grid size
    Real dx, dy, dz;
    Real *base_buffer, *u2, *u3; // Intermediate buffers
    C2Decomp *c2d; // Decomposition object
    bool periodicBC[3];
    int pRow, pCol;
    Real *fft_inputX, *fft_inputY, *fft_inputZ; // FFT input buffer
    Real *fft_outputX, *fft_outputY, *fft_outputZ; // FFT input buffer

    void performFFT();

    void performIFFT();

    void solveEigenvalues();

public:
    poissonSolver(index_t Nx, index_t Ny, index_t Nz, Real dx, Real dy, Real dz, C2Decomp *c2d)
        : Nx(Nx), Ny(Ny), Nz(Nz), dx(dx), dy(dy), dz(dz), c2d(c2d) {
        // Allocate buffers
        c2d->allocY(u2);
        c2d->allocZ(u3);
        // Retrieve grid dimensions
        
        // Allocate FFT buffers
        fft_inputX = fftwr_malloc(sizeof(Real) * c2d->xSize[0] * c2d->xSize[1] * c2d->xSize[2]);
        fft_outputX = fftwr_malloc(sizeof(Real) * c2d->xSize[0] * c2d->xSize[1] * c2d->xSize[2]);
        fft_inputY = fftwr_malloc(sizeof(Real) * c2d->ySize[0] * c2d->ySize[1] * c2d->ySize[2]);
        fft_outputY = fftwr_malloc(sizeof(Real) * c2d->ySize[0] * c2d->ySize[1] * c2d->ySize[2]);
        fft_inputZ = fftwr_malloc(sizeof(Real) * c2d->zSize[0] * c2d->zSize[1] * c2d->zSize[2]);
        fft_outputZ = fftwr_malloc(sizeof(Real) * c2d->zSize[0] * c2d->zSize[1] * c2d->zSize[2]);
    }

    ~poissonSolver() {
        fftwr_free(u2);
        fftwr_free(u3);
    }


    void setBuffer(Real *b); // Setter function to update `b`, the rhs
    void solve(Real *X);

    // later on change this to Real *X
    // Updates the input array `X` with the solution
};

#endif
