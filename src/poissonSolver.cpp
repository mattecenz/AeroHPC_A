#include "poissonSolver.hpp"

void poissonSolver::performFFT()
{
    // Perform FFT along x (for each y and z slice)
    for (int k = 0; k < c2d->xSize[2]; k++){
        for (int j = 0; j < c2d->xSize[1]; j++){
            const int base_index = c2d->xSize[0] * c2d->xSize[1] * k + c2d->xSize[0] * j;
            memcpy(fft_inputX, &base_buffer[base_index], sizeof(Real) * c2d->xSize[0]);
            fftwr_execute(fftwr_plan_r2r_1d(c2d->xSize[0], fft_inputX, fft_outputX, FFTW_REDFT10, FFTW_ESTIMATE)); // Execute FFT on each slice
            memcpy(&base_buffer[base_index], fft_outputX, sizeof(Real) * c2d->xSize[0]);
        }
    }

    c2d->transposeX2Y_MajorIndex(base_buffer, u2);

    for (int i = 0; i < c2d->ySize[0]; i++){
        for (int k = 0; k < c2d->ySize[2]; k++){
            const int base_index = c2d->ySize[2] * c2d->ySize[1] * i + c2d->ySize[1] * k;
            memcpy(fft_inputY, &u2[base_index], sizeof(Real) * c2d->ySize[1]);
            fftwr_execute(fftwr_plan_r2r_1d(c2d->ySize[1], fft_inputY, fft_outputY, FFTW_REDFT10, FFTW_ESTIMATE)); // Execute FFT on each slice
            memcpy(&u2[base_index], fft_outputY, sizeof(Real) * c2d->ySize[1]);
        }
    }

    c2d->transposeY2Z_MajorIndex(u2, u3);

    for (int j = 0; j < c2d->zSize[1]; j++){
        for (int i = 0; i < c2d->zSize[0]; i++){
            const int base_index = c2d->zSize[2] * c2d->zSize[0] * j + c2d->zSize[2] * i;
            memcpy(fft_inputZ, &u3[base_index], sizeof(Real) * c2d->zSize[2]);
            fftwr_execute(fftwr_plan_r2r_1d(c2d->zSize[2], fft_inputZ, fft_outputZ, FFTW_REDFT10, FFTW_ESTIMATE)); // Execute FFT on each slice
            memcpy(&u3[base_index], fft_outputZ, sizeof(Real) * c2d->zSize[2]);
        }
    }
}

void poissonSolver::performIFFT()
{

    for (int j = 0; j < c2d->zSize[1]; j++) {
        for (int i = 0; i < c2d->zSize[0]; i++) {
            const int base_index = c2d->zSize[2]*c2d->zSize[0] * j + c2d->zSize[2] * i;
            memcpy(fft_inputZ, &u3[base_index], sizeof(Real) * c2d->zSize[2]);
            fftwr_execute(fftwr_plan_r2r_1d(c2d->zSize[2], fft_inputZ, fft_outputZ, FFTW_REDFT01, FFTW_ESTIMATE));
            memcpy(&u3[base_index], fft_outputZ, sizeof(Real) * c2d->zSize[2]);
        }
    }

    c2d->transposeZ2Y_MajorIndex(u3, u2);


    for (int i = 0; i < c2d->ySize[0]; i++) {
        for (int k = 0; k < c2d->ySize[2]; k++) {
            const int base_index = c2d->ySize[2]*c2d->ySize[1] * i + c2d->ySize[1] * k;
            memcpy(fft_inputY, &u2[base_index], sizeof(Real) * c2d->ySize[1]);
            fftwr_execute(fftwr_plan_r2r_1d(c2d->ySize[1], fft_inputY, fft_outputY, FFTW_REDFT01, FFTW_ESTIMATE));
            memcpy(&u2[base_index], fft_outputY, sizeof(Real) * c2d->ySize[1]);
        }
    }

    c2d->transposeY2X_MajorIndex(u2, base_buffer);

    for (int k = 0; k < c2d->xSize[2]; k++) {
        for (int j = 0; j < c2d->xSize[1]; j++) {
            const int base_index = c2d->xSize[0]*c2d->xSize[1] * k + c2d->xSize[0] * j;
            memcpy(fft_inputX, &base_buffer[base_index], sizeof(Real) * c2d->xSize[0]);
            fftwr_execute(fftwr_plan_r2r_1d(c2d->xSize[0], fft_inputX, fft_outputX, FFTW_REDFT01, FFTW_ESTIMATE));
            memcpy(&base_buffer[base_index], fft_outputX, sizeof(Real) * c2d->xSize[0]);
        }
    }
}

void poissonSolver::solveEigenvalues()
{
    // Compute and apply eigenvalues
    for (int kp = 0; kp < c2d->zSize[2]; ++kp)
    {
        int kz = ((c2d->zStart[2] + kp) < Nz / 2) ? (c2d->zStart[2] + kp) : (c2d->zStart[2] + kp - Nz);
        for (int jp = 0; jp < c2d->zSize[1]; ++jp)
        {
            int ky = ((c2d->zStart[1] + jp) < Ny / 2) ? (c2d->zStart[1] + jp) : (c2d->zStart[1] + jp - Ny);
            for (int ip = 0; ip < c2d->zSize[0]; ++ip)
            {
                int kx = ((c2d->zStart[0] + ip) < Nx / 2) ? (c2d->zStart[0] + ip) : (c2d->zStart[0] + ip - Nx);
                int idx = kp * (c2d->zSize[1] * c2d->zSize[0]) + jp * c2d->zSize[0] + ip;
                Real eigenvalue = -(kx * kx + ky * ky + kz * kz);
                // im not sure whether i have to skip the zero eigen value or
                // divide by a very small number to simmulate zero
                if (eigenvalue != 0)
                {
                    u3[idx] /= eigenvalue;
                }
                else
                {
                    u3[idx] = 0.0; // Set zero mode to 0
                }
            }
        }
    }
}
void poissonSolver::solve(Real *buffer)
{
    // Initialize u1 with b (rhs)
    this->base_buffer = buffer; // Dynamically update `b`

    // perform forward FFT
    performFFT();
    // Solve eigenvalues
    solveEigenvalues();
    // Perform inverse FFT
    performIFFT();

    int scaling_fact = 8 * (Nx * Ny * Nz);
    for (int i = 0; i < c2d->xSize[0] * c2d->xSize[1] * c2d->xSize[2]; ++i)
    {
        // double check the normalization factor
        base_buffer[i] /= scaling_fact;
    }
}