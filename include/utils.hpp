/**
@file utils.h
@brief Math utilities
*/

#ifndef AEROHPC_A_UTILS_H
#define AEROHPC_A_UTILS_H

#include <StaggeredGrid.hpp>
namespace utils
{
    constexpr int CHANGE_THIS = 3;

    /**
    @param[in] grid Staggered grid
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the x direction
    */
    template <typename T, Addressing_T A>
    inline T d_dx(const StaggeredGrid<T, A> &grid, Component c, int i, int j, int k);

    /**
    @param[in] grid Staggered grid
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the y direction
    */
    template <typename T, Addressing_T A>
    inline T d_dy(const StaggeredGrid<T, A> &grid, Component c, int i, int j, int k);

    /**
    @param[in] grid Staggered grid
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return First derivative evaluation in point (i,j,k)
    @brief Computes the value of the first order derivative in a point along the z direction
    */
    template <typename T, Addressing_T A>
    inline T d_dz(const StaggeredGrid<T, A> &grid, Component c, int i, int j, int k);

    /**
    @param[in] grid Staggered grid
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the x direction
    */
    template <typename T, Addressing_T A>
    inline T d2_dx2(const StaggeredGrid<T, A> &grid, Component c, int i, int j, int k);

    /**
    @param[in] grid Staggered grid
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the y direction
    */
    template <typename T, Addressing_T A>
    inline T d2_dy2(const StaggeredGrid<T, A> &grid, Component c, int i, int j, int k);

    /**
    @param[in] grid Staggered grid
    @param[in] c Component whose derivative has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Second derivative evaluation in point (i,j,k)
    @brief Computes the value of the second order derivative in a point along the z direction
    */
    template <typename T, Addressing_T A>
    inline T d2_dz2(const StaggeredGrid<T, A> &grid, Component c, int i, int j, int k);


    /**
    @param[in] grid Staggered grid
    @param[in] c Component whose laplacian has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Laplacian evaluation in point (i,j,k)
    @brief Computes the value of the laplacian in a point
    */
    template <typename T, Addressing_T A>
    inline T lap(const StaggeredGrid<T, A> &grid, Component c, int i, int j, int k);

    /**
    @param[in] grid Staggered grid
    @param[in] to Destination grid of the interpolation
    @param[in] from Source grid of the interpolation
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Interpolated value
    @brief Computes the interpolation of a component from the grid "from" to the grid "to" in the point i,j,k of the "to" grid
    */
    template <typename T, Addressing_T A>
    inline T get_interpolation(const StaggeredGrid<T, A> &grid, Component to, Component from, int i, int j, int k);

    /**
    @param[in] grid Staggered grid
    @param[in] c Component whose convective term has to be computed
    @param[in] i,j,k Coordinates in the staggered grid corresponding to the component
    @return Convective term evaluation in point (i,j,k)
    @brief Computes the value of convective term of the Navier-Stokes equation in a point
    */
    template <typename T, Addressing_T A>
    inline T conv(const StaggeredGrid<T, A> &grid, Component c, int i, int j, int k);
}
#endif