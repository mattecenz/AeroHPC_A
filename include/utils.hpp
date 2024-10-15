#ifndef AEROHPC_A_UTILS_H
#define AEROHPC_A_UTILS_H

#include <StaggeredGrid.hpp>
namespace utils
{
    constexpr int CHANGE_THIS = 3;

    // Utils for first order derivatives
    template <typename T, Addressing_T A>
    inline T d_dx(const StaggeredGrid<T, A> &grid, Component c, int i, int j, int k)
    {
        return (grid(c, i + 1, j, k) - grid(c, i - 1, j, k)) / (CHANGE_THIS);
    }

    template <typename T, Addressing_T A>
    inline T d_dy(const StaggeredGrid<T, A> &grid, int i, int j, int k)
    {
        return (grid(c, i, j + 1, k) - grid(c, i, j - 1, k)) / (CHANGE_THIS);
    }

    template <typename T, Addressing_T A>
    inline T d_dz(const StaggeredGrid<T, A> &grid, int i, int j, int k)
    {
        return (grid(c, i, j, k + 1) - grid(c, i, j, k - 1)) / (CHANGE_THIS);
    }

    // Utils for second order derivatives
    template <typename T, Addressing_T A>
    inline T d2_dx2(const StaggeredGrid<T, A> &grid,Component c, int i, int j, int k)
    {
        return (grid(c, i - 1, j, k) - 2 * grid(c, i, j, k) + grid(c, i + 1, j, k)) / (CHANGE_THIS * CHANGE_THIS);
    }

    template <typename T, Addressing_T A>
    inline T d2_dy2(const StaggeredGrid<T, A> &grid, int i, int j, int k)
    {
        return (grid(c, i, j - 1, k) - 2 * grid(c, i, j, k) + grid(c, i, j + 1, k)) / (CHANGE_THIS * CHANGE_THIS);
    }

    template <typename T, Addressing_T A>
    inline T d2_dz2(const StaggeredGrid<T, A> &grid, int i, int j, int k)
    {
        return (grid(c, i, j, k - 1) - 2 * grid(c, i, j, k) + grid(c, i, j, k + 1)) / (CHANGE_THIS * CHANGE_THIS);
    }





    // Util for interpolation between staggered grids
    template <typename T, Addressing_T A>
    inline T get_interpolation(const StaggeredGrid<T, A> &grid, Component to, Component from, int i, int j, int k)
    {
        switch (to)
        {
        case Component::U:
            if (from == Component::V)
            {
                // Interpolate from V grid to U grid
                return (grid(Component::V, i, j, k) + grid(Component::V, i + 1, j, k) + grid(Component::V, i, j - 1, k) + grid(Component::V, i + 1, j - 1, k)) / 4;
            }
            else
            {
                // Interpolate from W grid to U grid
                return (grid(Component::W, i, j, k) + grid(Component::W, i + 1, j, k) + grid(Component::W, i, j, k - 1) + grid(Component::W, i + 1, j, k - 1)) / 4;
            }
            break;
        case Component::V:
            if (from == Component::U)
            {
                // Interpolate from U grid to V grid
                return (grid(Component::U, i, j, k) + grid(Component::U, i, j + 1, k) + grid(Component::U, i - 1, j, k) + grid(Component::U, i - 1, j + 1, k)) / 4;
            }
            else
            {
                // Interpolate from W grid to V grid
                return (grid(Component::W, i, j, k) + grid(Component::W, i, j + 1, k) + grid(Component::W, i, j, k - 1) + grid(Component::W, i, j + 1, k - 1)) / 4;
            }
            break;
        case Component::W:
            if (from == Component::U)
            {
                // Interpolate from U grid to W grid
                return (grid(Component::U, i, j, k) + grid(Component::U, i, j, k + 1) + grid(Component::U, i - 1, j, k) + grid(Component::U, i - 1, j, k + 1)) / 4;
            }
            else
            {
                // Interpolate from V grid to W grid
                return (grid(Component::V, i, j, k) + grid(Component::V, i, j, k + 1) + grid(Component::V, i, j - 1, k) + grid(Component::V, i, j - 1, k + 1)) / 4;
            }
            break;
        }
    }
}
#endif