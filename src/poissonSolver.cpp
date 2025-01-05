#include "poissonSolver.hpp"

poissonSolver::poissonSolver(int N, Real L, C2Decomp *c2d)
    : N(N), L(L), c2d(c2d) {
    dx = L / N;

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

    // Allocate FFT buffers
    fft_inputX  = (Real*)fftwr_malloc(sizeof(Real) * xSize[0] * xSize[1] * xSize[2]);
    fft_outputX = (Real*)fftwr_malloc(sizeof(Real) * xSize[0] * xSize[1] * xSize[2]);

    fft_inputY  = (Real*)fftwr_malloc(sizeof(Real) * ySize[0] * ySize[1] * ySize[2]);
    fft_outputY = (Real*)fftwr_malloc(sizeof(Real) * ySize[0] * ySize[1] * ySize[2]);

    fft_inputZ  = (Real*)fftwr_malloc(sizeof(Real) * zSize[0] * zSize[1] * zSize[2]);
    fft_outputZ = (Real*)fftwr_malloc(sizeof(Real) * zSize[0] * zSize[1] * zSize[2]);

}

poissonSolver::~poissonSolver() {
    fftwr_free(u1);
    fftwr_free(u2);
    fftwr_free(u3);
}

void poissonSolver::setB(Real *b) {
    this->b = b; // Dynamically update `b`
}


void poissonSolver::performFFT() {
    // Perform FFT along each axis

    // Perform FFT along x (for each y and z slice)
    for (int jp = 0; jp < xSize[1]; ++jp) {
        for (int kp = 0; kp < xSize[2]; ++kp) {
            for (int ip = 0; ip < xSize[0]; ++ip) {
                // int ii = ip * xSize[2] * xSize[1] + jp * xSize[2] + kp;
                int ii = kp * xSize[1] * xSize[0] + jp * xSize[0] + ip;

                fft_inputX[ip] = u1[ii];
            }
            fftwr_execute(fftwr_plan_r2r_1d(xSize[0], fft_inputX, fft_outputX, FFTW_REDFT10, FFTW_ESTIMATE));  // Execute FFT on each slice

            // Store the FFT result back into u1
            for (int ip = 0; ip < xSize[0]; ++ip) {
                // int ii = ip * xSize[2] * xSize[1] + jp * xSize[2] + kp;
                int ii = kp * xSize[1] * xSize[0] + jp * xSize[0] + ip;

                u1[ii] = fft_outputX[ip];
            }
        }
    }

    // Transpose u1 to u2 (y-direction)
    c2d->transposeX2Y(u1, u2);

    // Perform FFT along y (for each x and z slice)
    for (int kp = 0; kp < ySize[2]; ++kp) {
        for (int ip = 0; ip < ySize[0]; ++ip) {
            for (int jp = 0; jp < ySize[1]; ++jp) {
                // int ii = ip * ySize[2] * ySize[1] + jp * ySize[2] + kp;
                int ii = kp * ySize[1] * ySize[0] + jp * ySize[0] + ip;

                fft_inputY[jp] = u2[ii];
            }
            fftwr_execute(fftwr_plan_r2r_1d(ySize[0], fft_inputY, fft_outputY, FFTW_REDFT10, FFTW_ESTIMATE));  // Execute FFT on each slice

            // Store the FFT result back into u2
            for (int jp = 0; jp < ySize[1]; ++jp) {
                // int ii = ip * ySize[2] * ySize[1] + jp * ySize[2] + kp;
                int ii = kp * ySize[1] * ySize[0] + jp * ySize[0] + ip;

                u2[ii] = fft_outputY[jp];
            }
        }
    }

    // Transpose u2 to u3 (z-direction)
    c2d->transposeY2Z(u2, u3);

    // Perform FFT along z (for each x and y slice)
    for (int jp = 0; jp < zSize[1]; ++jp) {
        for (int ip = 0; ip < zSize[0]; ++ip) {
            for (int kp = 0; kp < zSize[2]; ++kp) {
                // int ii = ip * zSize[2] * zSize[1] + jp * zSize[2] + kp;
                int ii = kp * zSize[1] * zSize[0] + jp * zSize[0] + ip;

                fft_inputZ[kp] = u3[ii];
            }
            fftwr_execute(fftwr_plan_r2r_1d(zSize[0], fft_inputZ, fft_outputZ, FFTW_REDFT10, FFTW_ESTIMATE));  // Execute FFT on each slice

            // Store the FFT result back into u3
            for (int kp = 0; kp < zSize[2]; ++kp) {
                // int ii = ip * zSize[2] * zSize[1] + jp * zSize[2] + kp;
                int ii = kp * zSize[1] * zSize[0] + jp * zSize[0] + ip;

                u3[ii] = fft_outputZ[kp];
            }
        }
    }
}

void poissonSolver::performIFFT() {
    // Perform inverse FFT along each axis

    // Perform inverse FFT along Z direction (on u3)
    for (int jp = 0; jp < zSize[1]; ++jp) {
        for (int ip = 0; ip < zSize[0]; ++ip) {
            for (int kp = 0; kp < zSize[2]; ++kp) {
                // int ii = ip * zSize[2] * zSize[1] + jp * zSize[2] + kp;
                int ii = kp * zSize[1] * zSize[0] + jp * zSize[0] + ip;

                fft_inputZ[kp] = u3[ii];
            }
            // Execute inverse FFT
            fftwr_execute(fftwr_plan_r2r_1d(zSize[0], fft_inputZ, fft_outputZ, FFTW_REDFT01, FFTW_ESTIMATE));

            // Store the inverse FFT result back into u3
            for (int kp = 0; kp < zSize[2]; ++kp) {
                // int ii = ip * zSize[2] * zSize[1] + jp * zSize[2] + kp;
                int ii = kp * zSize[1] * zSize[0] + jp * zSize[0] + ip;

                u3[ii] = fft_outputZ[kp];
            }
        }
    }

    c2d->transposeZ2Y(u3, u2);

    // Perform inverse FFT along Y direction (on u2)
    for (int kp = 0; kp < ySize[2]; ++kp) {
        for (int ip = 0; ip < ySize[0]; ++ip) {
            for (int jp = 0; jp < ySize[1]; ++jp) {
                // int ii = ip * ySize[2] * ySize[1] + jp * ySize[2] + kp;
                int ii = kp * ySize[1] * ySize[0] + jp * ySize[0] + ip;

                fft_inputY[jp] = u2[ii];
            }
            // Execute inverse FFT
            fftwr_execute(fftwr_plan_r2r_1d(ySize[0], fft_inputY, fft_outputY, FFTW_REDFT01, FFTW_ESTIMATE));

            // Store the inverse FFT result back into u2
            for (int jp = 0; jp < ySize[1]; ++jp) {
                // int ii = ip * ySize[2] * ySize[1] + jp * ySize[2] + kp;
                int ii = kp * ySize[1] * ySize[0] + jp * ySize[0] + ip;

                u2[ii] = fft_outputY[jp];
            }
        }
    }

    c2d->transposeY2X(u2, u1);

    // Perform inverse FFT along X direction (on u1)
    for (int jp = 0; jp < xSize[1]; ++jp) {
        for (int kp = 0; kp < xSize[2]; ++kp) {
            for (int ip = 0; ip < xSize[0]; ++ip) {
                // int ii = ip * xSize[2] * xSize[1] + jp * xSize[2] + kp;
                int ii = kp * xSize[1] * xSize[0] + jp * xSize[0] + ip;

                fft_inputX[ip] = u1[ii];
            }
            // Execute inverse FFT
            fftwr_execute(fftwr_plan_r2r_1d(xSize[0], fft_inputX, fft_outputX, FFTW_REDFT01, FFTW_ESTIMATE));

            // Store the inverse FFT result back into u1
            for (int ip = 0; ip < xSize[0]; ++ip) {
                // int ii = ip * xSize[2] * xSize[1] + jp * xSize[2] + kp;
                int ii = kp * xSize[1] * xSize[0] + jp * xSize[0] + ip;

                u1[ii] = fft_outputX[ip];
            }
        }
    }

}

void poissonSolver::solveEigenvalues() {
    // Compute and apply eigenvalues
    for (int kp = 0; kp < zSize[2]; ++kp) {
        int kz = ((c2d->zStart[2] + kp) < N / 2) ? (c2d->zStart[2] + kp) : (c2d->zStart[2] + kp - N);
        for (int jp = 0; jp < zSize[1]; ++jp) {
            int ky = ((c2d->zStart[1] + jp) < N / 2) ? (c2d->zStart[1] + jp) : (c2d->zStart[1] + jp - N);
            for (int ip = 0; ip < zSize[0]; ++ip) {
                int kx = ((c2d->zStart[0] + ip) < N / 2) ? (c2d->zStart[0] + ip) : (c2d->zStart[0] + ip - N);

                int idx = kp * (zSize[1] * zSize[0]) + jp * zSize[0] + ip;
                Real eigenvalue = -(kx * kx + ky * ky + kz * kz);

                // im not sure whether i have to skip the zero eigen value or
                // divide by a very small number to simmulate zero
                if (eigenvalue != 0) {
                    u3[idx] /= eigenvalue;
                } else {
                    u3[idx] = 0.0; // Set zero mode to 0
                }

            }
        }
    }
}

void poissonSolver::solve(Real *X) {

    // Initialize u1 with b (rhs)
    for (int i = 0; i < xSize[0] * xSize[1] * xSize[2]; ++i){
        u1[i] = b[i];
    }

    // perform forward FFT
    performFFT();

    // Solve eigenvalues
    solveEigenvalues();

    // Perform inverse FFT
    performIFFT();

    // Store the result in X
    for (int i = 0; i < xSize[0] * xSize[1] * xSize[2]; ++i) {
        // double check the normalization factor
        int scaling_fact = 8*(N*N*N);
        X[i] = u1[i]/scaling_fact;
    }
}