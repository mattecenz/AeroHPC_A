#ifndef FFTDATA_HPP
#define FFTDATA_HPP

#include "Traits.hpp"
#include "data/Parameters.hpp"
#include "C2Decomp.hpp"
#include "fftw3.h"

// Rename the fftw library function based on Real type:
// "fftw_" functions are for double type
// "fftwf_" functions are for float type
// The renamed format is "fftwr_" (where r stands for Real)
#if REAL_USE_FLOAT
#define fftwr_malloc(a) real_p(fftwf_malloc(a))
#define fftwr_free fftwf_free
#define fftwr_plan fftwf_plan
#define fftwr_destroy_plan fftwf_destroy_plan
#define fftwr_plan_r2r_1d fftwf_plan_r2r_1d
#define fftwr_execute fftwf_execute
#else
#define fftwr_malloc(a) real_p(fftw_malloc(a))
#define fftwr_free fftw_free
#define fftwr_plan fftw_plan
#define fftwr_destroy_plan fftw_destroy_plan
#define fftwr_plan_r2r_1d fftw_plan_r2r_1d
#define fftwr_execute fftw_execute
#endif

#define FFT_NEUMANN FFTW_REDFT10
#define FFT_PERIODIC FFTW_R2HC
#define FFT_KIND(params, Axis)  params.periodic##Axis ? FFT_PERIODIC : FFT_I_NEUMANN

#define FFT_I_NEUMANN FFTW_REDFT01
#define FFT_I_PERIODIC FFTW_HC2R
#define FFT_I_KIND(params, Axis)  params.periodic##Axis ? FFT_I_PERIODIC : FFT_I_NEUMANN

class FFTData {
public:
    Real *inputX, *inputY, *inputZ; // FFT input buffer
    Real *outputX, *outputY, *outputZ; // FFT input buffer
    fftwr_plan planx, plany, planz, planx_i, plany_i, planz_i;

    Real *base_buffer, *u2, *u3, *eigs; // Intermediate buffers

    Real scalingFactor;

    FFTData(const Parameters &params, const C2Decomp &c2D, Real *base_buffer) {
        // Allocate FFT buffers
        inputX = fftwr_malloc(sizeof(Real) * c2D.xSize[0]);
        outputX = fftwr_malloc(sizeof(Real) * c2D.xSize[0]);
        inputY = fftwr_malloc(sizeof(Real) * c2D.ySize[1]);
        outputY = fftwr_malloc(sizeof(Real) * c2D.ySize[1]);
        inputZ = fftwr_malloc(sizeof(Real) * c2D.zSize[2]);
        outputZ = fftwr_malloc(sizeof(Real) * c2D.zSize[2]);

        planx = fftwr_plan_r2r_1d(c2D.xSize[0], inputX, outputX, FFT_KIND(params, X), FFTW_ESTIMATE);
        plany = fftwr_plan_r2r_1d(c2D.ySize[1], inputY, outputY, FFT_KIND(params, Y), FFTW_ESTIMATE);
        planz = fftwr_plan_r2r_1d(c2D.zSize[2], inputZ, outputZ, FFT_KIND(params, Z), FFTW_ESTIMATE);

        planz_i = fftwr_plan_r2r_1d(c2D.zSize[2], inputZ, outputZ, FFT_I_KIND(params, Z), FFTW_ESTIMATE);
        plany_i = fftwr_plan_r2r_1d(c2D.ySize[1], inputY, outputY, FFT_I_KIND(params, Y), FFTW_ESTIMATE);
        planx_i = fftwr_plan_r2r_1d(c2D.xSize[0], inputX, outputX, FFT_I_KIND(params, X), FFTW_ESTIMATE);

        this->base_buffer = base_buffer;
        c2D.allocY(u2);
        c2D.allocZ(u3);
        c2D.allocZ(eigs);

#define scal(params, Axis) (params.periodic##Axis ? 2 * params.glob_n##Axis : 1)

        scalingFactor = real(scal(params, X) * scal(params, Y) * scal(params, Z));

        computeEigs(params, c2D);
    }

    ~FFTData() {
        fftwr_free(u2);
        fftwr_free(u3);
        fftwr_free(eigs);

        fftwr_free(inputX);
        fftwr_free(inputY);
        fftwr_free(inputZ);
        fftwr_free(outputX);
        fftwr_free(outputY);
        fftwr_free(outputZ);

        fftwr_destroy_plan(planx);
        fftwr_destroy_plan(plany);
        fftwr_destroy_plan(planz);
        fftwr_destroy_plan(planx_i);
        fftwr_destroy_plan(plany_i);
        fftwr_destroy_plan(planz_i);
    }


    void computeEigs(const Parameters &params, const C2Decomp &c2D) const {
        std::cout << "-------------------------- EIGENVALUES -------------------------" << std::endl;
        int zDirXSize = c2D.zSize[2];
        int zDirYSize = c2D.zSize[0];
        int zDirZSize = c2D.zSize[1];
        int zDirXStart = c2D.zStart[2];
        int zDirYStart = c2D.zStart[0];
        int zDirZStart = c2D.zStart[1];

        // In Z transpose form we have Z on X-axis, Y on Z-axis and X on Y-axis

#define eig(i, delta, N) ( -pow( 2.0 / delta * std::sin( (real(i) * M_PI) / (2.0 * real(N)) ) , 2) )

        for (int k = 0; k < zDirZSize; k++) {
            const index_t layer_idx = k * zDirXSize * zDirYSize;
            const Real lambda_1 = eig(k + zDirZStart, params.dY, params.glob_nY);
            for (int j = 0; j < zDirYSize; j++) {
                const index_t row_idx = layer_idx + j * zDirXSize;
                const Real lambda_2 = eig(j + zDirYStart, params.dX, params.glob_nX);
                for (int i = 0; i < zDirXSize; i++) {
                    const Real lambda_3 = eig(i, params.dZ, params.glob_nZ);
                    Real e = lambda_1 + lambda_2 + lambda_3;
                    eigs[row_idx + i] = e;
                }
            }
        }
        std::cout << "--------------------------------------------------------------" << std::endl;
    }
};

#endif //FFTDATA_HPP
