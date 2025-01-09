#ifndef FFTDATA_HPP
#define FFTDATA_HPP

#include "Traits.hpp"
#include "Parameters.hpp"
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

class FFTData {
public:
    Real *inputX, *inputY, *inputZ; // FFT input buffer
    Real *outputX, *outputY, *outputZ; // FFT input buffer
    fftwr_plan planx, plany, planz, planx_i, plany_i, planz_i;

    Real *base_buffer, *u2, *u3, *eigs; // Intermediate buffers

    FFTData(const Parameters &params, const C2Decomp &c2D, Real *base_buffer) {
        // Allocate FFT buffers
        inputX = fftwr_malloc(sizeof(Real) * c2D.xSize[0] * c2D.xSize[1] * c2D.xSize[2]);
        outputX = fftwr_malloc(sizeof(Real) * c2D.xSize[0] * c2D.xSize[1] * c2D.xSize[2]);
        inputY = fftwr_malloc(sizeof(Real) * c2D.ySize[0] * c2D.ySize[1] * c2D.ySize[2]);
        outputY = fftwr_malloc(sizeof(Real) * c2D.ySize[0] * c2D.ySize[1] * c2D.ySize[2]);
        inputZ = fftwr_malloc(sizeof(Real) * c2D.zSize[0] * c2D.zSize[1] * c2D.zSize[2]);
        outputZ = fftwr_malloc(sizeof(Real) * c2D.zSize[0] * c2D.zSize[1] * c2D.zSize[2]);

        planx = fftwr_plan_r2r_1d(c2D.xSize[0], inputX, outputX, FFTW_REDFT10, FFTW_ESTIMATE);
        plany = fftwr_plan_r2r_1d(c2D.ySize[1], inputY, outputY, FFTW_REDFT10, FFTW_ESTIMATE);
        planz = fftwr_plan_r2r_1d(c2D.zSize[2], inputZ, outputZ, FFTW_REDFT10, FFTW_ESTIMATE);

        planz_i = fftwr_plan_r2r_1d(c2D.zSize[2], inputZ, outputZ, FFTW_REDFT01, FFTW_ESTIMATE);
        plany_i = fftwr_plan_r2r_1d(c2D.ySize[1], inputY, outputY, FFTW_REDFT01, FFTW_ESTIMATE);
        planx_i = fftwr_plan_r2r_1d(c2D.xSize[0], inputX, outputX, FFTW_REDFT01, FFTW_ESTIMATE);

        this->base_buffer = base_buffer;
        c2D.allocY(u2);
        c2D.allocZ(u3);
        c2D.allocZ(eigs);

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

#define eig(i, delta, N) (2.0 / delta * std::cos((i * M_PI) / (2.0 * N)))

    void computeEigs(const Parameters &params, const C2Decomp &c2D) const {
        for (int j = 0; j < c2D.zSize[1]; j++){
            const Real lambda_2 = eig(j + c2D.zStart[1], params.dY, params.glob_nY);
            const index_t base_index_1 = j * c2D.zSize[0] * c2D.zSize[2];

            for (int i = 0; i < c2D.zSize[0]; i++){
                const index_t base_index_2 = base_index_1 + i * c2D.zSize[2];
                const Real lambda_1 = eig(i + c2D.zStart[0], params.dX, params.glob_nX);

                for (int k = 0; k < c2D.zSize[2]; k++){
                    const Real lambda_3 = eig(k, params.dZ, params.glob_nZ);
                    eigs[base_index_2 + k] = (lambda_1 * lambda_1 + lambda_2 * lambda_2 + lambda_3 * lambda_3);
                }
            }
        }
    }
};

#endif //FFTDATA_HPP
