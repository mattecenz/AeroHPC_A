#ifndef SOLVERDATA_HPP
#define SOLVERDATA_HPP

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

    fftData_ptr = new FFTData(params, c2D, rkData.pressure_buffer_data);

    // TODO init BCs

}

// INDEXING MACRO
#define indexing(i,j,k) (i + j * params.loc_nX + k * params.loc_nX * params.loc_nY)
#define ghosted_indexing(i,j,k) (i + j * params.loc_gnX + k * params.loc_gnX * params.loc_gnY)

// COMPONENTS FOR GHOSTED BUFFERS
#define get_U_ghosted(data_ptr) (data)
#define get_V_ghosted(data_ptr) (&(data[params.grid_gndim]))
#define get_W_ghosted(data_ptr) (&(data[params.grid_gndim * 2]))

// COMPONENTS FOR BUFFERS
#define get_U(data) (data)
#define get_V(data) (&data[params.grid_ndim])
#define get_W(data) (&data[params.grid_ndim * 2])

// ELEMENTS FOR GHOSTED BUFFERS
#define U(data_ptr, i, j, k) get_U_ghosted(data_ptr)[ghosted_indexing(i,j,k)]
#define V(data_ptr, i, j, k) get_V_ghosted(data_ptr)[ghosted_indexing(i,j,k)]
#define W(data_ptr, i, j, k) get_W_ghosted(data_ptr)[ghosted_indexing(i,j,k)]
#define P(data_ptr, i, j, k) data_ptr[ghosted_indexing(i,j,k)]

// ELEMENTS FOR BUFFERS
#define pressure_buffer(i,j,k) pressure_buffer[indexing(i, j, k)]
#define rhs_U(i,j,k) get_U(rhs_buffer)[indexing(i, j, k)]
#define rhs_V(i,j,k) get_V(rhs_buffer)[indexing(i, j, k)]
#define rhs_W(i,j,k) get_W(rhs_buffer)[indexing(i, j, k)]

#endif //SOLVERDATA_HPP
