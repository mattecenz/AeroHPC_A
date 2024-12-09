#include "poissonSolver.hpp"

poissonSolver::poissonSolver(int N, double L, double *b, bool periodicBC[3], int pRow, int pCol)
    : N(N), L(L), b(b), pRow(pRow), pCol(pCol) {
    dx = L / N;

    // Initialize periodic boundary conditions
    this->periodicBC[0] = periodicBC[0];
    this->periodicBC[1] = periodicBC[1];
    this->periodicBC[2] = periodicBC[2];

    // Initialize decomposition
    c2d = new C2Decomp(N, N, N, pRow, pCol, this->periodicBC);

    // Allocate buffers
    c2d->allocX(u1);
    c2d->allocY(u2);
    c2d->allocZ(u3);

    // Retrieve grid dimensions
    xSize[0] = c2d->xSize[0];
    xSize[1] = c2d->xSize[1];
    xSize[2] = c2d->xSize[2];
    ySize[0] = c2d->ySize[0];
    ySize[1] = c2d->ySize[1];
    ySize[2] = c2d->ySize[2];
    zSize[0] = c2d->zSize[0];
    zSize[1] = c2d->zSize[1];
    zSize[2] = c2d->zSize[2];
}

poissonSolver::~poissonSolver() {
    delete c2d;

    fftw_destroy_plan(fft_plan_x);
    fftw_destroy_plan(fft_plan_y);
    fftw_destroy_plan(fft_plan_z);
    fftw_destroy_plan(ifft_plan_x);
    fftw_destroy_plan(ifft_plan_y);
    fftw_destroy_plan(ifft_plan_z);

    fftw_free(u1);
    fftw_free(u2);
    fftw_free(u3);
}

void poissonSolver::initializeGrid() {
    // Distribute b into u1
    for (int kp = 0; kp < xSize[2]; kp++) {
        for (int jp = 0; jp < xSize[1]; jp++) {
            for (int ip = 0; ip < xSize[0]; ip++) {
                int ii = kp * xSize[1] * xSize[0] + jp * xSize[0] + ip;

                int global_index = (c2d->xStart[2] + kp) * (N * N) 
                                 + (c2d->xStart[1] + jp) * N 
                                 + (c2d->xStart[0] + ip);

                u1[ii] = b[global_index];
            }
        }
    }
}

void poissonSolver::performFFT() {
    // Perform FFT along each axis
    c2d->fftX(u1, u1);
    c2d->transposeX2Y(u1, u2);
    c2d->fftY(u2, u2);
    c2d->transposeY2Z(u2, u3);
    c2d->fftZ(u3, u3);
}

void poissonSolver::solveEigenvalues() {
    // Compute and apply eigenvalues
    for (int kp = 0; kp < zSize[2]; ++kp) {
        int kz = ((c2d->zStart[2] + kp) < N / 2) ? (c2d->zStart[2] + kp) : (c2d->zStart[2] + kp - N);
        for (int jp = 0; jp < zSize[1]; ++jp) {
            int ky = ((c2d->zStart[1] + jp) < N / 2) ? (c2d->zStart[1] + jp) : (c2d->zStart[1] + jp - N);
            for (int ip = 0; ip < zSize[0]; ++ip) {
                int kx = ((c2d->zStart[0] + ip) < N / 2) ? (c2d->zStart[0] + ip) : (c2d->zStart[0] + ip - N);

                size_t idx = kp * zSize[1] * zSize[0] + jp * zSize[0] + ip;
                double ksq = -(kx * kx + ky * ky + kz * kz);
                double eigenvalue = (ksq == 0) ? 1e-10 : ksq;
                u3[idx] /= eigenvalue;
            }
        }
    }
}

void poissonSolver::performIFFT() {
    // Perform inverse FFT along each axis
    c2d->ifftZ(u3, u3);
    c2d->transposeZ2Y(u3, u2);
    c2d->ifftY(u2, u2);
    c2d->transposeY2X(u2, u1);
    c2d->ifftX(u1, u1);
}

void poissonSolver::solve(double *X) {
    initializeGrid();
    performFFT();
    solveEigenvalues();
    performIFFT();

    // Store the result in X
    for (int i = 0; i < xSize[0] * xSize[1] * xSize[2]; ++i) {
        X[i] = u1[i];
    }
}
