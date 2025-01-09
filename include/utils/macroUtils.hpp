#ifndef MACROUTILS_HPP
#define MACROUTILS_HPP

#define ITERATE_DOMAIN_VELOCITY(i,j,k)                                      \
const index_t nz = params.isOnRight ? params.loc_nZ - 1 : params.loc_nZ;    \
const index_t ny = params.isOnTop ? params.loc_nY - 1 : params.loc_nY;      \
const index_t nx = params.loc_nY - 1;                                       \
for (index_t k = 0; k < nz; ++k) {                                          \
    for (index_t j = 0; j < ny; ++j) {                                      \
    _Pragma("omp simd")                                                     \
        for (index_t i = 0; i < nx; ++i) {

#define ITERATE_DOMAIN_PRESSURE(i,j,k)                                  \
const index_t nz = params.loc_nZ;                                       \
const index_t ny = params.loc_nY;                                       \
const index_t nx = params.loc_nX;                                       \
for (index_t k = 0; k < nz; ++k) {                                      \
    for (index_t j = 0; j < ny; ++j) {                                  \
    _Pragma("omp simd")                                                 \
        for (index_t i = 0; i < nx; ++i) {

#define ITERATE_DOMAIN_END()                                        \
        }                                                           \
    }                                                               \
}

#endif //MACROUTILS_HPP
