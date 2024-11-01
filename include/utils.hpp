/**
@file utils.h
@brief Math utilities
*/

#ifndef AEROHPC_A_UTILS_H
#define AEROHPC_A_UTILS_H

#include "Grid.hpp"

namespace utils {
    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the x direction
    */
#define d_dx(C) inline Real d_dx_##C(const Grid &grid, int i, int j, int k){ \
    return (grid.C(i + 1, j, k) - grid.C(i - 1, j, k)) / (2 * grid.dx); \
}

    d_dx(U)

    d_dx(V)

    d_dx(W)

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the y direction
    */
#define d_dy(C) inline Real d_dy_##C(const Grid &grid, int i, int j, int k){ \
    return (grid.C(i, j + 1, k) - grid.C(i, j - 1, k)) / (2 * grid.dy); \
}

    d_dy(U)

    d_dy(V)

    d_dy(W)

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the z direction
    */
#define d_dz(C) inline Real d_dz_##C(const Grid &grid, int i, int j, int k){ \
    return (grid.C(i, j, k + 1) - grid.C(i, j, k - 1)) / (2 * grid.dz); \
}

    d_dz(U)

    d_dz(V)

    d_dz(W)

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the x direction
    */
#define d2_dx2(C) inline Real d2_dx2_##C(const Grid &grid, int i, int j, int k){ \
    return (grid.C(i - 1, j, k) - 2 * grid.C(i, j, k) + grid.C(i + 1, j, k)) / (grid.dx * grid.dx); \
}

    d2_dx2(U)

    d2_dx2(V)

    d2_dx2(W)

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the y direction
    */
#define d2_dy2(C) inline Real d2_dy2_##C(const Grid &grid, int i, int j, int k){ \
    return (grid.C(i, j - 1, k) - 2 * grid.C(i, j, k) + grid.C(i, j + 1, k)) / (grid.dy * grid.dy); \
}

    d2_dy2(U)

    d2_dy2(V)

    d2_dy2(W)

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the z direction
    */
#define d2_dz2(C) inline Real d2_dz2_##C(const Grid &grid, int i, int j, int k){ \
    return (grid.C(i, j, k - 1) - 2 * grid.C(i, j, k) + grid.C(i, j, k + 1)) / (grid.dz * grid.dz); \
}

    d2_dz2(U)

    d2_dz2(V)

    d2_dz2(W)


    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose laplacian has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Laplacian evaluation in point (i,j,k)
    @brief Computes the value of the laplacian in a point
    */
#define lap(C) inline Real lap_##C(const Grid &grid, int i, int j, int k){ \
    return d2_dx2_##C(grid, i, j, k) + d2_dy2_##C(grid, i, j, k) + d2_dz2_##C(grid, i, j, k); \
}

    lap(U)

    lap(V)

    lap(W)

    /**
    @param[in] grid Staggered grid
    @param[in] to Destination grid of the interpolation
    @param[in] from Source grid of the interpolation
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Interpolated value
    @brief Computes the interpolation of a component from the grid "from" to the grid "to" in the point i,j,k of the "to" grid
    */
#define interp(to, from) inline Real intp_##from##_on_##to(const Grid &grid, const int i, const int j, const int k) { \
    const int U = 0; \
    const int V = 1; \
    const int W = 2; \
    if constexpr (to == 0) { \
        if constexpr (from == 1) { \
            return (grid.V(i, j, k) + grid.V(i + 1, j, k) + \
                    grid.V(i, j - 1, k) + grid.V(i + 1, j - 1, k)) / 4; \
        } else { \
            return (grid.W(i, j, k) + grid.W(i + 1, j, k) + \
                    grid.W(i, j, k - 1) + grid.W(i + 1, j, k - 1)) / 4; \
        } \
    } else if constexpr (to == 1) { \
        if constexpr (from == 0) { \
            return (grid.U(i, j, k) + grid.U(i, j + 1, k) + \
                    grid.U(i - 1, j, k) + grid.U(i - 1, j + 1, k)) / 4; \
        } else { \
            return (grid.W(i, j, k) + grid.W(i, j + 1, k) + \
                    grid.W(i, j, k - 1) + grid.W(i, j + 1, k - 1)) / 4; \
        } \
    } else { \
        if constexpr (from == 0) { \
            return (grid.U(i, j, k) + grid.U(i, j, k + 1) + \
                    grid.U(i - 1, j, k) + grid.U(i - 1, j, k + 1)) / 4; \
        } else { \
            return (grid.V(i, j, k) + grid.V(i, j, k + 1) + \
                    grid.V(i, j - 1, k) + grid.V(i, j - 1, k + 1)) / 4; \
        } \
    } \
}

    interp(U, V)

    interp(U, W)

    interp(V, U)

    interp(V, W)

    interp(W, U)

    interp(W, V)

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose convective term has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Convective term evaluation in point (i,j,k)
    @brief Computes the value of convective term of the Navier-Stokes equation in a point
    */
#define conv(C) inline Real conv_##C(const Grid &grid, const int i, const int j, const int k) { \
    const int U = 0; \
    const int V = 1; \
    const int W = 2; \
    if constexpr (C == 0) { \
        return grid.U(i, j, k) * d_dx_U(grid, i, j, k) + \
               intp_V_on_U(grid, i, j, k) * d_dy_U(grid, i, j, k) + \
               intp_W_on_U(grid, i, j, k) * d_dz_U(grid, i, j, k); \
    } else if constexpr (C == 1) { \
        return intp_U_on_V(grid, i, j, k) * d_dx_V(grid, i, j, k) + \
               grid.V(i, j, k) * d_dy_V(grid, i, j, k) + \
               intp_W_on_V(grid, i, j, k) * d_dz_V(grid, i, j, k); \
    } else { \
        return intp_U_on_W(grid, i, j, k) * d_dx_W(grid, i, j, k) + \
               intp_V_on_W(grid, i, j, k) * d_dy_W(grid, i, j, k) + \
               grid.W(i, j, k) * d_dz_W(grid, i, j, k); \
    } \
}

    conv(U)

    conv(V)

    conv(W)

}
#endif