#include <utils.hpp>
using namespace utils;

namespace utils
{

    template <>
    Real d_dx(const StaggeredGrid<Real, Addressing_T::STANDARD> &grid, Component c, int i, int j, int k)
    {
        return (grid(c, i + 1, j, k) - grid(c, i - 1, j, k)) / (CHANGE_THIS);
    }

    template <>
    Real d_dy(const StaggeredGrid<Real, Addressing_T::STANDARD> &grid, Component c, int i, int j, int k)
    {
        return (grid(c, i, j + 1, k) - grid(c, i, j - 1, k)) / (CHANGE_THIS);
    }

    template <>
    Real d_dz(const StaggeredGrid<Real, Addressing_T::STANDARD> &grid, Component c, int i, int j, int k)
    {
        return (grid(c, i, j, k + 1) - grid(c, i, j, k - 1)) / (CHANGE_THIS);
    }

    template <>
    Real d2_dx2(const StaggeredGrid<Real, Addressing_T::STANDARD> &grid, Component c, int i, int j, int k)
    {
        return (grid(c, i - 1, j, k) - 2 * grid(c, i, j, k) + grid(c, i + 1, j, k)) / (CHANGE_THIS * CHANGE_THIS);
    }

    template <>
    Real d2_dy2(const StaggeredGrid<Real, Addressing_T::STANDARD> &grid, Component c, int i, int j, int k)
    {
        return (grid(c, i, j - 1, k) - 2 * grid(c, i, j, k) + grid(c, i, j + 1, k)) / (CHANGE_THIS * CHANGE_THIS);
    }

    template <>
    Real d2_dz2(const StaggeredGrid<Real, Addressing_T::STANDARD> &grid, Component c, int i, int j, int k)
    {
        return (grid(c, i, j, k - 1) - 2 * grid(c, i, j, k) + grid(c, i, j, k + 1)) / (CHANGE_THIS * CHANGE_THIS);
    }

    template <>
    Real get_interpolation(const StaggeredGrid<Real, Addressing_T::STANDARD> &grid, Component to, Component from, int i, int j, int k)
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