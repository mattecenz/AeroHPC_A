#ifndef PRESSURESOLVER_HPP
#define PRESSURESOLVER_HPP

#include "SolverData.hpp"

#define xDirXSize c2D.xSize[0]
#define xDirYSize c2D.xSize[1]
#define xDirZSize c2D.xSize[2]
#define yDirXSize c2D.ySize[0]
#define yDirYSize c2D.ySize[1]
#define yDirZSize c2D.ySize[2]
#define zDirXSize c2D.zSize[0]
#define zDirYSize c2D.zSize[1]
#define zDirZSize c2D.zSize[2]

inline void computeFFT(){
    // Perform FFT along x (for each y and z slice)
    for (int k = 0; k < xDirYSize; k++){
        for (int j = 0; j < xDirYSize; j++){
            const int base_index = xDirXSize * xDirYSize * k + xDirYSize * j;
            memcpy(fftData.inputX, &fftData.base_buffer[base_index], sizeof(Real) * xDirXSize);
            fftwr_execute(fftData.planx); // Execute FFT on each slice
            memcpy(&fftData.base_buffer[base_index], fftData.outputX, sizeof(Real) * xDirXSize);
        }
    }

    c2D.transposeX2Y_MajorIndex(fftData.base_buffer, fftData.u2);

    for (int i = 0; i < yDirXSize; i++){
        for (int k = 0; k < yDirZSize; k++){
            const int base_index = yDirXSize * yDir * i + c2d->ySize[1] * k;
            memcpy(fft_inputY, &u2[base_index], sizeof(Real) * c2d->ySize[1]);
            fftwr_execute(plany); // Execute FFT on each slice
            memcpy(&u2[base_index], fft_outputY, sizezof(Real) * c2d->ySize[1]);
        }
    }

    c2d->transposeY2Z_MajorIndex(u2, u3);

    for (int j = 0; j < c2d->zSize[1]; j++){
        for (int i = 0; i < c2d->zSize[0]; i++){
            const int base_index = c2d->zSize[2] * c2d->zSize[0] * j + c2d->zSize[2] * i;
            memcpy(fft_inputZ, &u3[base_index], sizeof(Real) * c2d->zSize[2]);
            fftwr_execute(planz); // Execute FFT on each slice
            memcpy(&u3[base_index], fft_outputZ, sizeof(Real) * c2d->zSize[2]);
        }
    }
}


#endif //PRESSURESOLVER_HPP
