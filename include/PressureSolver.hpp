#ifndef PRESSURESOLVER_HPP
#define PRESSURESOLVER_HPP

#include "data/SolverData.hpp"

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
    for (int k = 0; k < xDirZSize; k++){
        for (int j = 0; j < xDirYSize; j++){
            const int row_idx = xDirXSize * xDirYSize * k + xDirXSize * j;
            memcpy(fftData.inputX, &fftData.base_buffer[row_idx], sizeof(Real) * xDirXSize);
            fftwr_execute(fftData.planx); // Execute FFT on each slice
            memcpy(&fftData.base_buffer[row_idx], fftData.outputX, sizeof(Real) * xDirXSize);
        }
    }

    c2D.transposeX2Y_MajorIndex(fftData.base_buffer, fftData.u2);

    for (int i = 0; i < yDirXSize; i++){
        for (int k = 0; k < yDirZSize; k++){
            const int row_idx = yDirZSize * yDirYSize * i + yDirYSize * k;
            memcpy(fftData.inputY, &fftData.u2[row_idx], sizeof(Real) * yDirYSize);
            fftwr_execute(fftData.plany); // Execute FFT on each slice
            memcpy(&fftData.u2[row_idx], fftData.outputY, sizeof(Real) * yDirYSize);
        }
    }

    c2D.transposeY2Z_MajorIndex(fftData.u2, fftData.u3);

    for (int j = 0; j < zDirYSize; j++){
        for (int i = 0; i < zDirXSize; i++){
            const int row_idx = zDirZSize * zDirXSize * j + zDirZSize * i;
            memcpy(fftData.inputZ, &fftData.u3[row_idx], sizeof(Real) * zDirZSize);
            fftwr_execute(fftData.planz); // Execute FFT on each slice
            memcpy(&fftData.u3[row_idx], fftData.outputZ, sizeof(Real) * zDirZSize);
        }
    }
}

inline void computeIFFT()
{
    for (int j = 0; j < zDirYSize; j++) {
        for (int i = 0; i < zDirXSize; i++) {
            const int row_idx = zDirZSize * zDirXSize * j + zDirZSize * i;
            memcpy(fftData.inputZ, &fftData.u3[row_idx], sizeof(Real) * zDirZSize);
            fftwr_execute(fftData.planz_i);
            memcpy(&fftData.u3[row_idx], fftData.outputZ, sizeof(Real) * zDirZSize);
        }
    }

    c2D.transposeZ2Y_MajorIndex(fftData.u3, fftData.u2);


    for (int i = 0; i < yDirXSize; i++) {
        for (int k = 0; k < yDirZSize; k++) {
            const int row_idx = yDirYSize * yDirZSize * i + yDirYSize * k;
            memcpy(fftData.inputY, &fftData.u2[row_idx], sizeof(Real) * yDirYSize);
            fftwr_execute(fftData.plany_i);
            memcpy(&fftData.u2[row_idx], fftData.outputY, sizeof(Real) * yDirYSize);
        }
    }

    c2D.transposeY2X_MajorIndex(fftData.u2, fftData.base_buffer);

    for (int k = 0; k < xDirZSize; k++) {
        for (int j = 0; j < xDirYSize; j++) {
            const int row_idx = xDirXSize * xDirYSize * k + xDirXSize * j;
            memcpy(fftData.inputX, &fftData.base_buffer[row_idx], sizeof(Real) * xDirXSize);
            fftwr_execute(fftData.planx_i);
            memcpy(&fftData.base_buffer[row_idx], fftData.outputX, sizeof(Real) * xDirXSize);
        }
    }
}

inline void solveEigenvalues(){
    for (int j = 0; j < zDirYSize; j++){
        const index_t layer_idx = j * zDirXSize * zDirZSize;
        for (int i = 0; i < zDirXSize; i++){
            const index_t row_idx = layer_idx + i * zDirZSize;
            for (int k = 0; k < zDirZSize; k++){
                if(fftData.eigs[row_idx+k]!=0.0){
                    fftData.u3[row_idx + k] /= fftData.eigs[row_idx + k];
                }else{
                    fftData.u3[row_idx + k] = 0.0;
                }
            }
        }
    }
}

inline void solvePressure(){
    computeFFT();
    solveEigenvalues();
    computeIFFT();

    for(index_t idx = 0; idx < xDirXSize * xDirYSize * xDirZSize; ++idx){
        fftData.base_buffer[idx] /= fftData.scalingFactor;
    }
}




#endif //PRESSURESOLVER_HPP
