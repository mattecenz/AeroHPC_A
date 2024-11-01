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
     Real d_dx(const Grid<A> &grid, Component c, int i, int j, int k);

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the y direction
    */
    template<Addressing_T A>
     Real d_dy(const Grid<A> &grid, Component c, int i, int j, int k);

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the z direction
    */
    template<Addressing_T A>
     Real d_dz(const Grid<A> &grid, Component c, int i, int j, int k);

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the x direction
    */
    template<Addressing_T A>
     Real d2_dx2(const Grid<A> &grid, Component c, int i, int j, int k);

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the y direction
    */
    template<Addressing_T A>
     Real d2_dy2(const Grid<A> &grid, Component c, int i, int j, int k);

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the z direction
    */
    template<Addressing_T A>
     Real d2_dz2(const Grid<A> &grid, Component c, int i, int j, int k);


    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose laplacian has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Laplacian evaluation in point (i,j,k)
    @brief Computes the value of the laplacian in a point
    */
    template<Addressing_T A>
     Real lap(const Grid<A> &grid, Component c, int i, int j, int k);

    template<Addressing_T A>
     Vector lap(const Grid<A> &grid, int i, int j, int k);

    /**
    @param[in] grid Staggered grid
    @param[in] to Destination grid of the interpolation
    @param[in] from Source grid of the interpolation
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Interpolated value
    @brief Computes the interpolation of a component from the grid "from" to the grid "to" in the point i,j,k of the "to" grid
    */
    template<Addressing_T A>
     Real get_interpolation(const Grid<A> &grid, Component to, Component from, int i, int j, int k);

    /**
    @param[in] model Model containing both the staggered grid and the spacing
    @param[in] c Component whose convective term has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Convective term evaluation in point (i,j,k)
    @brief Computes the value of convective term of the Navier-Stokes equation in a point
    */
    template<Addressing_T A>
     Real conv(const Grid<A> &grid, Component c, int i, int j, int k);

    template<Addressing_T A>
     Vector conv(const Grid<A> &grid, int i, int j, int k);
}
#endif