#include "FFT.hpp"

inline void FFT::computeFFT()
{

}

inline void FFT::computeIFFT()
{

    for (int j = 0; j < c2d->zSize[1]; j++) {
        for (int i = 0; i < c2d->zSize[0]; i++) {
            const int base_index = c2d->zSize[2]*c2d->zSize[0] * j + c2d->zSize[2] * i;
            memcpy(fft_inputZ, &u3[base_index], sizeof(Real) * c2d->zSize[2]);
            fftwr_execute(planz_i);
            memcpy(&u3[base_index], fft_outputZ, sizeof(Real) * c2d->zSize[2]);
        }
    }

    c2d->transposeZ2Y_MajorIndex(u3, u2);


    for (int i = 0; i < c2d->ySize[0]; i++) {
        for (int k = 0; k < c2d->ySize[2]; k++) {
            const int base_index = c2d->ySize[2]*c2d->ySize[1] * i + c2d->ySize[1] * k;
            memcpy(fft_inputY, &u2[base_index], sizeof(Real) * c2d->ySize[1]);
            fftwr_execute(plany_i);
            memcpy(&u2[base_index], fft_outputY, sizeof(Real) * c2d->ySize[1]);
        }
    }

    c2d->transposeY2X_MajorIndex(u2, base_buffer);

    for (int k = 0; k < c2d->xSize[2]; k++) {
        for (int j = 0; j < c2d->xSize[1]; j++) {
            const int base_index = c2d->xSize[0]*c2d->xSize[1] * k + c2d->xSize[0] * j;
            memcpy(fft_inputX, &base_buffer[base_index], sizeof(Real) * c2d->xSize[0]);
            fftwr_execute(planx_i);
            memcpy(&base_buffer[base_index], fft_outputX, sizeof(Real) * c2d->xSize[0]);
        }
    }
}



inline void FFT::solveEigenvalues()
{
  	for (int j = 0; j < c2d->zSize[1]; j++){
        const index_t base_index_1 = j * c2d->zSize[0] * c2d->zSize[2];
        for (int i = 0; i < c2d->zSize[0]; i++){
            const index_t base_index_2 = base_index_1 + i * c2d->zSize[2];
            for (int k = 0; k < c2d->zSize[2]; k++){
                u3[base_index_2 + k] /= eigs[base_index_2 + k];
            }
        }
    }
}


void FFT::solve()
{
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

