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
    template<Addressing_T A>
    inline Real d_dx(const Grid<A> &grid, Component c, int i, int j, int k) {
        return (grid(c, i + 1, j, k) - grid(c, i - 1, j, k)) / (2 * grid.dx);
    }

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the y direction
    */
    template<Addressing_T A>
    inline Real d_dy(const Grid<A> &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j + 1, k) - grid(c, i, j - 1, k)) / (2 * grid.dy);
    }

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the z direction
    */
    template<Addressing_T A>
    inline Real d_dz(const Grid<A> &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j, k + 1) - grid(c, i, j, k - 1)) / (2 * grid.dz);
    }

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the x direction
    */
    template<Addressing_T A>
    inline Real d2_dx2(const Grid<A> &grid, Component c, int i, int j, int k) {
        return (grid(c, i - 1, j, k) - 2 * grid(c, i, j, k) + grid(c, i + 1, j, k)) /
               (grid.dx * grid.dx);
    }

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the y direction
    */
    template<Addressing_T A>
    inline Real d2_dy2(const Grid<A> &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j - 1, k) - 2 * grid(c, i, j, k) + grid(c, i, j + 1, k)) /
               (grid.dy * grid.dy);
    }

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the z direction
    */
    template<Addressing_T A>
    inline Real d2_dz2(const Grid<A> &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j, k - 1) - 2 * grid(c, i, j, k) + grid(c, i, j, k + 1)) /
               (grid.dz * grid.dz);
    }

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose laplacian has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Laplacian evaluation in point (i,j,k)
    @brief Computes the value of the laplacian in a point
    */
    template<Addressing_T A>
    inline Real lap(const Grid<A> &grid, Component c, int i, int j, int k) {
        return d2_dx2(grid, c, i, j, k) + d2_dy2(grid, c, i, j, k) + d2_dz2(grid, c, i, j, k);
    }

    template<Addressing_T A>
    inline Vector lap(const Grid<A> &grid, int i, int j, int k) {
        return {
            d2_dx2(grid, U, i, j, k) + d2_dy2(grid, U, i, j, k) + d2_dz2(grid, U, i, j, k),
            d2_dx2(grid, V, i, j, k) + d2_dy2(grid, V, i, j, k) + d2_dz2(grid, V, i, j, k),
            d2_dx2(grid, W, i, j, k) + d2_dy2(grid, W, i, j, k) + d2_dz2(grid, W, i, j, k)
        };
    }

    /**
    @param[in] grid Staggered grid
    @param[in] to Destination grid of the interpolation
    @param[in] from Source grid of the interpolation
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Interpolated value
    @brief Computes the interpolation of a component from the grid "from" to the grid "to" in the point i,j,k of the "to" grid
    */
    template<Addressing_T A>
    inline Real get_interpolation(const Grid<A> &grid, Component to, Component from, int i, int j, int k) {
        switch (to) {
            case Component::U:
                if (from == Component::V) {
                    return (grid(Component::V, i, j, k) + grid(Component::V, i + 1, j, k) +
                            grid(Component::V, i, j - 1, k) + grid(Component::V, i + 1, j - 1, k)) / 4;
                } else {
                    return (grid(Component::W, i, j, k) + grid(Component::W, i + 1, j, k) +
                            grid(Component::W, i, j, k - 1) + grid(Component::W, i + 1, j, k - 1)) / 4;
                }
            case Component::V:
                if (from == Component::U) {
                    return (grid(Component::U, i, j, k) + grid(Component::U, i, j + 1, k) +
                            grid(Component::U, i - 1, j, k) + grid(Component::U, i - 1, j + 1, k)) / 4;
                } else {
                    return (grid(Component::W, i, j, k) + grid(Component::W, i, j + 1, k) +
                            grid(Component::W, i, j, k - 1) + grid(Component::W, i, j + 1, k - 1)) / 4;
                }
            case Component::W:
                if (from == Component::U) {
                    return (grid(Component::U, i, j, k) + grid(Component::U, i, j, k + 1) +
                            grid(Component::U, i - 1, j, k) + grid(Component::U, i - 1, j, k + 1)) / 4;
                } else {
                    return (grid(Component::V, i, j, k) + grid(Component::V, i, j, k + 1) +
                            grid(Component::V, i, j - 1, k) + grid(Component::V, i, j - 1, k + 1)) / 4;
                }
            default:
                return 0;
        }
    }

    /**
    @param[in] grid Staggered grid
    @param[in] c Component whose convective term has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Convective term evaluation in point (i,j,k)
    @brief Computes the value of the convective term for component c in point (i,j,k)
    */
    template<Addressing_T A>
    inline Real conv(const Grid<A> &grid, Component c, int i, int j, int k) {
        switch (c) {
            case Component::U:
                return grid(c, i, j, k) * d_dx(grid, c, i, j, k) +
                       get_interpolation(grid, Component::U, Component::V, i, j, k) * d_dy(grid, c, i, j, k) +
                       get_interpolation(grid, Component::U, Component::W, i, j, k) * d_dz(grid, c, i, j, k);
            case Component::V:
                return get_interpolation(grid, Component::V, Component::U, i, j, k) * d_dx(grid, c, i, j, k) +
                       grid(c, i, j, k) * d_dy(grid, c, i, j, k) +
                       get_interpolation(grid, Component::V, Component::W, i, j, k) * d_dz(grid, c, i, j, k);
            case Component::W:
                return get_interpolation(grid, Component::W, Component::U, i, j, k) * d_dx(grid, c, i, j, k) +
                       get_interpolation(grid, Component::W, Component::V, i, j, k) * d_dy(grid, c, i, j, k) +
                       grid(c, i, j, k) * d_dz(grid, c, i, j, k);
            default:
                return 0;
        }
    }

    template<Addressing_T A>
    inline Vector conv(const Grid<A> &grid, int i, int j, int k) {
        return {
            grid(U, i, j, k) * d_dx(grid, U, i, j, k) +
                       get_interpolation(grid, U, V, i, j, k) * d_dy(grid, U, i, j, k) +
                       get_interpolation(grid, U, W, i, j, k) * d_dz(grid, U, i, j, k),
            get_interpolation(grid, V, U, i, j, k) * d_dx(grid, V, i, j, k) +
                       grid(V, i, j, k) * d_dy(grid, V, i, j, k) +
                       get_interpolation(grid, V, W, i, j, k) * d_dz(grid, V, i, j, k),
            get_interpolation(grid, W, U, i, j, k) * d_dx(grid, W, i, j, k) +
                       get_interpolation(grid, W, V, i, j, k) * d_dy(grid, W, i, j, k) +
                       grid(W, i, j, k) * d_dz(grid, W, i, j, k)
        };
    }
} // namespace utils

#endif // AEROHPC_A_UTILS_H
