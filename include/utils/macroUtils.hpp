#ifndef MACROUTILS_HPP
#define MACROUTILS_HPP

#define NO_SKIP 0
#define SKIP_U 1
#define SKIP_V 2
#define SKIP_W 4

#define ITERATE_DOMAIN_VELOCITY(i, j, k, compute_physics, skip_C)                           \
{                                                                                           \
    const bool skip_U = skip_C & SKIP_U;                                                    \
    const bool skip_V = skip_C & SKIP_V;                                                    \
    const bool skip_W = skip_C & SKIP_W;                                                    \
    const index_t nz = (params.isOnRight && skip_W) ? params.loc_nZ - 1 : params.loc_nZ;    \
    const index_t ny = (params.isOnTop && skip_V) ? params.loc_nY - 1 : params.loc_nY;      \
    const index_t nx = skip_U ? params.loc_nX - 1 : params.loc_nX;                          \
    for (index_t k = 0; k < nz; ++k) {                                                      \
        Real z = 0;                                                                         \
        if constexpr (compute_physics)                                                      \
            z = real(k + params.st_nZ) * params.dZ + params.originZ;                        \
        for (index_t j = 0; j < ny; ++j) {                                                  \
            Real y = 0;                                                                     \
            if constexpr (compute_physics)                                                  \
                y = real(j + params.st_nY) * params.dY + params.originY;                    \
        _Pragma("omp simd")                                                                 \
            for (index_t i = 0; i < nx; ++i) {                                              \
                Real x = 0;                                                                 \
                if constexpr (compute_physics)                                              \
                    x = real(i + params.st_nX) * params.dX + params.originX;

#define ITERATE_DOMAIN_PRESSURE(i, j, k, compute_physics)                               \
{                                                                                       \
    const index_t nz = params.loc_nZ;                                                   \
    const index_t ny = params.loc_nY;                                                   \
    const index_t nx = params.loc_nX;                                                   \
    for (index_t k = 0; k < nz; ++k) {                                                  \
        Real z = 0;                                                                     \
        if constexpr (compute_physics)                                                  \
            z = real(k + params.st_nZ) * params.dZ + params.originZ;                    \
        for (index_t j = 0; j < ny; ++j) {                                              \
            Real y = 0;                                                                 \
            if constexpr (compute_physics)                                              \
                y = real(j + params.st_nY) * params.dY + params.originY;                \
        _Pragma("omp simd")                                                             \
            for (index_t i = 0; i < nx; ++i) {                                          \
                Real x = 0;                                                             \
                if constexpr (compute_physics)                                          \
                    x = real(i + params.st_nX) * params.dX + params.originX;


#define ITERATE_DOMAIN_END()    \
            }                   \
        }                       \
    }                           \
}

#define ITERATE_X_ROW(i,x)                                                  \
    for (index_t i = 0; i < params.loc_nX; i++) {                           \
        const Real x = real(i + params.st_nX) * params.dX + params.originX;

#define ITERATE_XZ_FACE(i,k,x,z)                                                \
    for (index_t k = 0; k < params.loc_nZ; k++) {                               \
        const Real z = real(k + params.st_nZ) * params.dZ + params.originZ;     \
        ITERATE_X_ROW(i,x)

#define ITERATE_XY_FACE(i,j,x,y)                                                \
    for (index_t j = 0; j < params.loc_nY; j++) {                               \
        const Real y = real(j + params.st_nY) * params.dY + params.originY;     \
        ITERATE_X_ROW(i,x)

#define ITERATE_YZ_FACE(j,k,y,z)                                                \
    for (index_t j = 0; j < params.loc_nY; j++) {                               \
        const Real y = real(j + params.st_nY) * params.dY + params.originY;     \
        for (index_t k = 0; k < params.loc_nZ; k++) {                           \
            const Real z = real(k + params.st_nZ) * params.dZ + params.originZ;

#define ITERATE_FACE_END()  \
    }                       \
}


#endif //MACROUTILS_HPP
