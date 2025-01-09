#ifndef SOLVERINFO_HPP
#define SOLVERINFO_HPP

#include "Parameters.hpp"
#include "C2Decomp.hpp"
#include "RKData.hpp"
#include "FFTData.hpp"

#define params_ptr _params
#define c2D_ptr _c2D
#define rkData_ptr _rkData
#define fftData_ptr _fftData

#define params (*params_ptr)
#define c2D (*c2D_ptr)
#define rkData (*rkData_ptr)
#define fftData (*fftData_ptr)

inline Parameters params = nullptr;
inline C2Decomp c2D = nullptr;
inline RKData rkData = nullptr;
inline FFTData fftData = nullptr;

void inline destroyData() {
    delete params_ptr;
    delete c2D_ptr;
    delete rkData_ptr;
    delete fftData_ptr;
}

void inline initData(const Real dimX, const Real dimY, const Real dimZ,
                     const Real originX, const Real originY, const Real originZ,
                     const Real deltaT, const index_t timestep, const Real Re,
                     const index_t nX, const index_t nY, const index_t nZ,
                     const index_t nPY, const index_t nPZ,
                     const bool periodicBC[3]) {

    c2D_ptr = new C2Decomp(nX, nY, nZ, nPY, nPZ, periodicBC);

    params_ptr = new Parameters(dimX, dimY, dimZ, originX, originY, originZ,
                             nX, nY, nZ, deltaT, timestep, Re, c2D);

    rkData_ptr = new RKData(params);

    fftData_ptr = new FFTData(params, c2D, data.pressure_buffer);

    // TODO init BCs

}

#define indexing(i,j,k) i + j * params.loc_nX + k * params.loc_nX * params.loc_nY
#define ghosted_indexing(i,j,k) i + j * params.loc_gnX + k * params.loc_gnX * params.loc_gnY


// TODO MODIFY
#define velocity(i,j,k) velocity_data[ghosted_indexing(i, j, k)]
#define pressure(i,j,k) pressure_data[ghosted_indexing(i, j, k)]

#define velocity_buffer(i,j,k) velocity_buffer[ghosted_indexing(i, j, k)]
#define pressure_buffer(i,j,k) pressure_buffer[indexing(i, j, k)]
#define rhs_buffer(i,j,k) rhs_buffer[indexing(i, j, k)]

#endif //SOLVERDATA_HPP
