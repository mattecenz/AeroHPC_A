#ifndef PRESSURESOLVER_HPP
#define PRESSURESOLVER_HPP

#include "data/SolverData.hpp"

#define xDirXSize c2D.xSize[0]
#define xDirYSize c2D.xSize[1]
#define xDirZSize c2D.xSize[2]
#define yDirXSize c2D.ySize[1]
#define yDirYSize c2D.ySize[2]
#define yDirZSize c2D.ySize[0]
#define zDirXSize c2D.zSize[2]
#define zDirYSize c2D.zSize[0]
#define zDirZSize c2D.zSize[1]

inline void computeFFT() {
    //**********/ FFT FORWARD ON X AXIS /******************************************************************//
    for (int k = 0; k < xDirZSize; k++) {
        for (int j = 0; j < xDirYSize; j++) {
            const int row_idx = xDirXSize * xDirYSize * k + xDirXSize * j;
            memcpy(fftData.inputX, &fftData.base_buffer[row_idx], sizeof(Real) * xDirXSize);
            fftwr_execute(fftData.planx); // Execute FFT on each slice
            memcpy(&fftData.base_buffer[row_idx], fftData.outputX, sizeof(Real) * xDirXSize);
        }
    }
    //*****************************************************************************************************//

    // GET Y REPRESENTATION
    c2D.transposeX2Y_MajorIndex(fftData.base_buffer, fftData.u2);

    //**********/ FFT FORWARD ON Y AXIS /******************************************************************//
    for (int k = 0; k < yDirZSize; k++) {
        for (int j = 0; j < yDirYSize; j++) {
            const int row_idx = yDirXSize * yDirYSize * k + yDirXSize * j;
            memcpy(fftData.inputY, &fftData.u2[row_idx], sizeof(Real) * yDirXSize);
            fftwr_execute(fftData.plany); // Execute FFT on each slice
            memcpy(&fftData.u2[row_idx], fftData.outputY, sizeof(Real) * yDirXSize);
        }
    }
    //*****************************************************************************************************//

    // GET Z REPRESENTATION
    c2D.transposeY2Z_MajorIndex(fftData.u2, fftData.u3);

    //**********/ FFT FORWARD ON Z AXIS /******************************************************************//
    for (int k = 0; k < zDirZSize; k++) {
        for (int j = 0; j < zDirYSize; j++) {
            const int row_idx = zDirXSize * zDirYSize * k + zDirXSize * j;
            memcpy(fftData.inputZ, &fftData.u3[row_idx], sizeof(Real) * zDirXSize);
            fftwr_execute(fftData.planz); // Execute FFT on each slice
            memcpy(&fftData.u3[row_idx], fftData.outputZ, sizeof(Real) * zDirXSize);
        }
    }
    //*****************************************************************************************************//
}

inline void computeIFFT() {
    //**********/ FFT BACKWARD ON Z AXIS /*****************************************************************//
    for (int k = 0; k < zDirZSize; k++) {
        for (int j = 0; j < zDirYSize; j++) {
            const int row_idx = zDirXSize * zDirYSize * k + zDirXSize * j;
            memcpy(fftData.inputZ, &fftData.u3[row_idx], sizeof(Real) * zDirXSize);
            fftwr_execute(fftData.planz_i);
            memcpy(&fftData.u3[row_idx], fftData.outputZ, sizeof(Real) * zDirXSize);
        }
    }
    //*****************************************************************************************************//

    // GET Y REPRESENTATION
    c2D.transposeZ2Y_MajorIndex(fftData.u3, fftData.u2);

    //**********/ FFT BACKWARD ON Y AXIS /*****************************************************************//
    for (int k = 0; k < yDirZSize; k++) {
        for (int j = 0; j < yDirYSize; j++) {
            const int row_idx = yDirXSize * yDirYSize * k + yDirXSize * j;
            memcpy(fftData.inputY, &fftData.u2[row_idx], sizeof(Real) * yDirXSize);
            fftwr_execute(fftData.plany_i);
            memcpy(&fftData.u2[row_idx], fftData.outputY, sizeof(Real) * yDirXSize);
        }
    }
    //*****************************************************************************************************//

    // GET X REPRESENTATION
    c2D.transposeY2X_MajorIndex(fftData.u2, fftData.base_buffer);

    //**********/ FFT BACKWARD ON X AXIS /*****************************************************************//
    for (int k = 0; k < xDirZSize; k++) {
        for (int j = 0; j < xDirYSize; j++) {
            const int row_idx = xDirXSize * xDirYSize * k + xDirXSize * j;
            memcpy(fftData.inputX, &fftData.base_buffer[row_idx], sizeof(Real) * xDirXSize);
            fftwr_execute(fftData.planx_i);
            memcpy(&fftData.base_buffer[row_idx], fftData.outputX, sizeof(Real) * xDirXSize);
        }
    }
    //*****************************************************************************************************//
}

inline void solveEigenvalues() {
    // APPLY EIGENVALUES TO Z REPRESENTATION
    for (int k = 0; k < zDirZSize; k++) {
        const index_t layer_idx = k * zDirXSize * zDirYSize;
        for (int j = 0; j < zDirYSize; j++) {
            const index_t row_idx = layer_idx + j * zDirXSize;
            for (int i = 0; i < zDirXSize; i++) {
                fftData.u3[row_idx + i] *= fftData.eigs[row_idx + i];
            }
        }
    }
}

inline void solvePressure() {
    computeFFT();
    solveEigenvalues();
    computeIFFT();
    /*
    // APPLY SCALING FACTOR TO X REPRESENTATION
    for (index_t idx = 0; idx < xDirXSize * xDirYSize * xDirZSize; ++idx) {
        fftData.base_buffer[idx] /= fftData.scalingFactor;
    }
    */
}


#endif //PRESSURESOLVER_HPP
