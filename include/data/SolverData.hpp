#ifndef SOLVERDATA_HPP
#define SOLVERDATA_HPP

#include "data/Parameters.hpp"
#include "data/RKData.hpp"
#include "data/FFTData.hpp"
#include "data/Constants.hpp"
#include "data/DomainData.hpp"
#include "C2Decomp.hpp"

#define params_ptr _params
#define c2D_ptr _c2D
#define rkData_ptr _rkData
#define fftData_ptr _fftData
#define consts_ptr _consts
#define domData_ptr _domData

#define params (*params_ptr)
#define c2D (*c2D_ptr)
#define rkData (*rkData_ptr)
#define fftData (*fftData_ptr)
#define consts (*consts_ptr)
#define domData (*domData_ptr)

inline Parameters params = nullptr;
inline C2Decomp c2D = nullptr;
inline RKData rkData = nullptr;
inline FFTData fftData = nullptr;
inline Constants consts = nullptr;
inline DomainData domData = nullptr;

// INDEXING MACRO
#define indexing(i,j,k) (i + j * params.loc_nX + k * params.loc_nX * params.loc_nY)
#define ghosted_indexing(i,j,k) ((i+1) + (j+1) * params.loc_gnX + (k+1) * params.loc_gnX * params.loc_gnY)

// COMPONENTS FOR GHOSTED BUFFERS
#define get_U_ghosted(data_ptr) (&data_ptr[0])
#define get_V_ghosted(data_ptr) (&(data_ptr[params.grid_gndim]))
#define get_W_ghosted(data_ptr) (&(data_ptr[params.grid_gndim * 2]))
#define get_P_ghosted(data_ptr) (&(data_ptr[params.grid_gndim * 3]))

// COMPONENTS FOR BUFFERS
#define get_U(data_ptr) (&data_ptr[0])
#define get_V(data_ptr) (&data_ptr[params.grid_ndim])
#define get_W(data_ptr) (&data_ptr[params.grid_ndim * 2])
#define get_P(data_ptr) (&data_ptr[params.grid_ndim * 3])

// ELEMENTS FOR GHOSTED BUFFERS
#define U(data_ptr, i, j, k) get_U_ghosted(data_ptr)[ghosted_indexing(i,j,k)]
#define V(data_ptr, i, j, k) get_V_ghosted(data_ptr)[ghosted_indexing(i,j,k)]
#define W(data_ptr, i, j, k) get_W_ghosted(data_ptr)[ghosted_indexing(i,j,k)]
#define P(data_ptr, i, j, k) get_P_ghosted(data_ptr)[ghosted_indexing(i,j,k)]

// ELEMENTS FOR BUFFERS
#define rhs_U(i,j,k) get_U(rkData.rhs_data)[indexing(i, j, k)]
#define rhs_V(i,j,k) get_V(rkData.rhs_data)[indexing(i, j, k)]
#define rhs_W(i,j,k) get_W(rkData.rhs_data)[indexing(i, j, k)]
#define rhs_P(i,j,k) get_P(rkData.rhs_data)[indexing(i, j, k)]

void inline initData(const Real dimX, const Real dimY, const Real dimZ,
                     const Real originX, const Real originY, const Real originZ,
                     const Real deltaT, const index_t timestep, const Real Re,
                     const index_t nX, const index_t nY, const index_t nZ,
                     const index_t nPY, const index_t nPZ,
                     const bool periodicBC[3],
                     DomainData *domainData) {

    c2D_ptr = new C2Decomp(nX, nY, nZ, nPY, nPZ, periodicBC);

    params_ptr = new Parameters(dimX, dimY, dimZ, originX, originY, originZ,
                             nX, nY, nZ, deltaT, timestep, Re, periodicBC, c2D);

    rkData_ptr = new RKData(params);

    fftData_ptr = new FFTData(params, c2D, get_P(rkData.rhs_data));

    consts_ptr = new Constants(Re, deltaT);

    _domData = domainData;
}

void inline destroyData() {
    delete params_ptr;
    delete c2D_ptr;
    delete rkData_ptr;
    delete fftData_ptr;
    delete consts_ptr;
    domData_ptr = nullptr;
}

#endif //SOLVERDATA_HPP
