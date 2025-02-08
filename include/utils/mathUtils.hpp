/**
@file mathUtils.hpp
@brief Math utilities
*/

#ifndef AEROHPC_A_MATHUTILS_H
#define AEROHPC_A_MATHUTILS_H

#include "data/SolverData.hpp"

namespace mathUtils {


#define interpolate_on_grid(C) \
    inline Real interp_##C##_onGrid(const Real *data, const index_t i, const index_t j, const index_t k) { \
        const index_t U = 0; \
        const index_t V = 1; \
        const index_t W = 2; \
        const index_t P = 3; \
        if constexpr (C == U) { \
            return (U(data, i - 1, j, k) + U(data, i - 1, j - 1, k) + U(data, i - 1, j, k - 1) + U(data, i - 1, j - 1, k - 1)) / 4; \
        } else if constexpr (C == V) {  \
            return (V(data, i, j - 1, k) + V(data, i - 1, j - 1, k) + V(data, i, j - 1, k - 1) + V(data, i - 1, j - 1, k - 1)) / 4; \
        } else if constexpr (C == W){ \
            return (W(data, i, j, k - 1) + W(data, i - 1, j, k - 1) + W(data, i, j - 1, k - 1) + W(data, i - 1, j - 1, k - 1)) / 4; \
        } else if constexpr (C == P){ \
            return (P(data, i, j, k) + P(data, i, j, k - 1) + P(data, i, j - 1, k) + P(data, i, j - 1, k - 1)  \
                    + P(data, i - 1, j, k) + P(data, i - 1, j, k - 1) + P(data, i - 1, j - 1, k) + P(data, i - 1, j - 1, k - 1)) / 8;  \
        } \
    }

    interpolate_on_grid(U)

    interpolate_on_grid(V)

    interpolate_on_grid(W)

    interpolate_on_grid(P)
#undef interpolate_on_grid


#define interpolate(C) inline Real intp_##C(const Real *data, const int i, const int j, const int k) { \
    const int U = 0; \
    const int V = 1; \
    const int W = 2; \
    if constexpr (C == U) { \
        return (U(data, i, j, k) + U(data, i - 1, j, k)) / 2; \
    } else if constexpr (C == V) { \
        return (V(data, i, j, k) + V(data, i, j - 1, k)) / 2; \
    } else { \
        return (W(data, i, j, k) + W(data, i, j, k - 1)) / 2; \
    } \
}

    interpolate(U)

    interpolate(V)

    interpolate(W)
#undef interpolate


    inline Real vel_div(const Real *data, int i, int j, int k) {
        return (U(data, i, j, k) - U(data, i - 1, j, k)) / params.dX
               + (V(data, i, j, k) - V(data, i, j - 1, k)) / params.dY
               + (W(data, i, j, k) - W(data, i, j, k - 1)) / params.dZ;
    }

    inline Real dp_dx_U(const Real *data, int i, int j, int k) {
        return (P(data, i + 1, j, k) - P(data, i, j, k)) / params.dX;
    }

    inline Real dp_dy_V(const Real *data, int i, int j, int k) {
        return (P(data, i, j + 1, k) - P(data, i, j, k)) / params.dY;
    }

    inline Real dp_dz_W(const Real *data, int i, int j, int k) {
        return (P(data, i, j, k + 1) - P(data, i, j, k)) / params.dZ;
    }


    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the x direction
    */
#define d_dx(C) inline Real d_dx_##C(const Real *data, int i, int j, int k){ \
    return (C(data, i + 1, j, k) - C(data, i - 1, j, k)) / (2 * params.dX); \
}

    d_dx(U)

    d_dx(V)

    d_dx(W)
#undef d_dx

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the y direction
    */
#define d_dy(C) inline Real d_dy_##C(const Real *data, int i, int j, int k){ \
    return (C(data, i, j + 1, k) - C(data, i, j - 1, k)) / (2 * params.dY); \
}

    d_dy(U)

    d_dy(V)

    d_dy(W)
#undef d_dy

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the z direction
    */
#define d_dz(C) inline Real d_dz_##C(const Real *data, int i, int j, int k){ \
    return (C(data, i, j, k + 1) - C(data, i, j, k - 1)) / (2 * params.dZ); \
}

    d_dz(U)

    d_dz(V)

    d_dz(W)
#undef d_dz

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the x direction
    */
#define d2_dx2(C) inline Real d2_dx2_##C(const Real *data, int i, int j, int k){ \
    return (C(data, i - 1, j, k) - 2 * C(data, i, j, k) + C(data, i + 1, j, k)) / (params.dX * params.dX); \
}

    d2_dx2(U)

    d2_dx2(V)

    d2_dx2(W)
#undef d2_dx2

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the y direction
    */
#define d2_dy2(C) inline Real d2_dy2_##C(const Real *data, int i, int j, int k){ \
    return (C(data, i, j - 1, k) - 2 * C(data, i, j, k) + C(data, i, j + 1, k)) / (params.dY * params.dY); \
}

    d2_dy2(U)

    d2_dy2(V)

    d2_dy2(W)
#undef d2_dy2

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the z direction
    */
#define d2_dz2(C) inline Real d2_dz2_##C(const Real *data, int i, int j, int k){ \
    return (C(data, i, j, k - 1) - 2 * C(data, i, j, k) + C(data, i, j, k + 1)) / (params.dZ * params.dZ); \
}

    d2_dz2(U)

    d2_dz2(V)

    d2_dz2(W)
#undef d2_dz2


    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose laplacian has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Laplacian evaluation in point (i,j,k)
    @brief Computes the value of the laplacian in a point
    */
#define lap(C) inline Real lap_##C(const Real *data, int i, int j, int k){ \
    return d2_dx2_##C(data, i, j, k) + d2_dy2_##C(data, i, j, k) + d2_dz2_##C(data, i, j, k); \
}

    lap(U)

    lap(V)

    lap(W)
#undef lap

    /**
    @param[in] grid Staggered grid
    @param[in] to Destination grid of the interpolation
    @param[in] from Source grid of the interpolation
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Interpolated value
    @brief Computes the interpolation of a component from the grid "from" to the grid "to" in the point i,j,k of the "to" grid
    */
#define interp(to, from) inline Real intp_##from##_on_##to(const Real *data, const int i, const int j, const int k) { \
    const int U = 0; \
    const int V = 1; \
    const int W = 2; \
    if constexpr (to == 0) { \
        if constexpr (from == 1) { \
            return (V(data,i, j, k) + V(data,i + 1, j, k) + \
                    V(data,i, j - 1, k) + V(data,i + 1, j - 1, k)) / 4; \
        } else { \
            return (W(data,i, j, k) + W(data,i + 1, j, k) + \
                    W(data,i, j, k - 1) + W(data,i + 1, j, k - 1)) / 4; \
        } \
    } else if constexpr (to == 1) { \
        if constexpr (from == 0) { \
            return (U(data,i, j, k) + U(data,i, j + 1, k) + \
                    U(data,i - 1, j, k) + U(data,i - 1, j + 1, k)) / 4; \
        } else { \
            return (W(data,i, j, k) + W(data,i, j + 1, k) + \
                    W(data,i, j, k - 1) + W(data,i, j + 1, k - 1)) / 4; \
        } \
    } else { \
        if constexpr (from == 0) { \
            return (U(data,i, j, k) + U(data,i, j, k + 1) + \
                    U(data,i - 1, j, k) + U(data,i - 1, j, k + 1)) / 4; \
        } else { \
            return (V(data,i, j, k) + V(data,i, j, k + 1) + \
                    V(data,i, j - 1, k) + V(data,i, j - 1, k + 1)) / 4; \
        } \
    } \
}

    interp(U, V)

    interp(U, W)

    interp(V, U)

    interp(V, W)

    interp(W, U)

    interp(W, V)
#undef interp

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose convective term has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Convective term evaluation in point (i,j,k)
    @brief Computes the value of convective term of the Navier-Stokes equation in a point
    */
#define conv(C) inline Real conv_##C(const Real *data, const int i, const int j, const int k) { \
    const int U = 0; \
    const int V = 1; \
    const int W = 2; \
    if constexpr (C == 0) { \
        return U(data, i, j, k) * d_dx_U(data, i, j, k) + \
               intp_V_on_U(data, i, j, k) * d_dy_U(data, i, j, k) + \
               intp_W_on_U(data, i, j, k) * d_dz_U(data, i, j, k); \
    } else if constexpr (C == 1) { \
        return intp_U_on_V(data, i, j, k) * d_dx_V(data, i, j, k) + \
               V(data, i, j, k) * d_dy_V(data, i, j, k) + \
               intp_W_on_V(data, i, j, k) * d_dz_V(data, i, j, k); \
    } else { \
        return intp_U_on_W(data, i, j, k) * d_dx_W(data, i, j, k) + \
               intp_V_on_W(data, i, j, k) * d_dy_W(data, i, j, k) + \
               W(data, i, j, k) * d_dz_W(data, i, j, k); \
    } \
}

    conv(U)

    conv(V)

    conv(W)
#undef conv

}
#endif //AEROHPC_A_MATHUTILS_H
