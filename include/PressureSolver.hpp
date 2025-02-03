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
    // Perform FFT along x (for each y and z slice)
    // std::cout << std::setprecision(3) << std::fixed;
    for (int k = 0; k < xDirZSize; k++) {
        for (int j = 0; j < xDirYSize; j++) {
            const int row_idx = xDirXSize * xDirYSize * k + xDirXSize * j;
            // for (int i = 0; i < xDirXSize; i++) {
            //     std::cout << fftData.base_buffer[row_idx + i] << " ";
            // }
            // std::cout << std::endl;
            memcpy(fftData.inputX, &fftData.base_buffer[row_idx], sizeof(Real) * xDirXSize);


            fftwr_execute(fftData.planx); // Execute FFT on each slice
            memcpy(&fftData.base_buffer[row_idx], fftData.outputX, sizeof(Real) * xDirXSize);
        }
        // std::cout << std::endl;
    }
    // std::cout << "=====================" << std::endl;


    c2D.transposeX2Y_MajorIndex(fftData.base_buffer, fftData.u2);


    for (int k = 0; k < yDirZSize; k++) {
        for (int j = 0; j < yDirYSize; j++) {
            const int row_idx = yDirXSize * yDirYSize * k + yDirXSize * j;
            // for (int i = 0; i < yDirXSize; i++) {
            //     std::cout << fftData.u2[row_idx + i] << " ";
            // }
            // std::cout << std::endl;
            memcpy(fftData.inputY, &fftData.u2[row_idx], sizeof(Real) * yDirXSize);
            fftwr_execute(fftData.plany); // Execute FFT on each slice
            memcpy(&fftData.u2[row_idx], fftData.outputY, sizeof(Real) * yDirXSize);
        }
        // std::cout << std::endl;
    }
    // std::cout << "=====================" << std::endl;


    c2D.transposeY2Z_MajorIndex(fftData.u2, fftData.u3);

    for (int k = 0; k < zDirZSize; k++) {
        for (int j = 0; j < zDirYSize; j++) {
            const int row_idx = zDirXSize * zDirYSize * k + zDirXSize * j;
            // for (int i = 0; i < zDirXSize; i++) {
            //     std::cout << fftData.u3[row_idx + i] << " ";
            // }
            // std::cout << std::endl;
            memcpy(fftData.inputZ, &fftData.u3[row_idx], sizeof(Real) * zDirXSize);
            fftwr_execute(fftData.planz); // Execute FFT on each slice
            memcpy(&fftData.u3[row_idx], fftData.outputZ, sizeof(Real) * zDirXSize);
        }
        // std::cout << std::endl;
    }
    // std::cout << "=====================" << std::endl;

}

inline void computeIFFT() {
    // std::cout << "========== INVERSO ===========" << std::endl;

    for (int k = 0; k < zDirZSize; k++) {
        for (int j = 0; j < zDirYSize; j++) {
            const int row_idx = zDirXSize * zDirYSize * k + zDirXSize * j;
            // for (int i = 0; i < zDirXSize; i++) {
            //     std::cout << fftData.u3[row_idx + i] << " ";
            // }
            // std::cout << std::endl;
            memcpy(fftData.inputZ, &fftData.u3[row_idx], sizeof(Real) * zDirXSize);
            fftwr_execute(fftData.planz_i);
            memcpy(&fftData.u3[row_idx], fftData.outputZ, sizeof(Real) * zDirXSize);
        }
        // std::cout << std::endl;
    }
    // std::cout << "=====================" << std::endl;

    c2D.transposeZ2Y_MajorIndex(fftData.u3, fftData.u2);

    for (int k = 0; k < yDirZSize; k++) {
        for (int j = 0; j < yDirYSize; j++) {
            const int row_idx = yDirXSize * yDirYSize * k + yDirXSize * j;
            // for (int i = 0; i < yDirXSize; i++) {
            //     std::cout << fftData.u2[row_idx + i] << " ";
            // }
            // std::cout << std::endl;
            memcpy(fftData.inputY, &fftData.u2[row_idx], sizeof(Real) * yDirXSize);
            fftwr_execute(fftData.plany_i);
            memcpy(&fftData.u2[row_idx], fftData.outputY, sizeof(Real) * yDirXSize);
        }
        // std::cout << std::endl;
    }
    // std::cout << "=====================" << std::endl;

    c2D.transposeY2X_MajorIndex(fftData.u2, fftData.base_buffer);

    for (int k = 0; k < xDirZSize; k++) {
        for (int j = 0; j < xDirYSize; j++) {
            const int row_idx = xDirXSize * xDirYSize * k + xDirXSize * j;
            // for (int i = 0; i < xDirXSize; i++) {
            //     std::cout << fftData.base_buffer[row_idx + i] << " ";
            // }
            // std::cout << std::endl;
            memcpy(fftData.inputX, &fftData.base_buffer[row_idx], sizeof(Real) * xDirXSize);
            fftwr_execute(fftData.planx_i);
            memcpy(&fftData.base_buffer[row_idx], fftData.outputX, sizeof(Real) * xDirXSize);
        }
        // std::cout << std::endl;
    }
}

inline void solveEigenvalues() {
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

    // std::cout << "------------------------Start of Output of the poisson------------------------" << std::endl;
    for (index_t idx = 0; idx < xDirXSize * xDirYSize * xDirZSize; ++idx) {
        // fftData.base_buffer[idx] *= fftData.scalingFactor;
        // std::cout << fftData.base_buffer[idx] << std::endl;
    }
    // std::cout << "-------------------------End of Output of the poisson-------------------------" << std::endl;
}


#endif //PRESSURESOLVER_HPP