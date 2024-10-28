#include <utils.hpp>

using namespace utils;

namespace utils {

    template<>
    Real d_dx(const Grid<STANDARD> &grid, Component c, int i, int j, int k) {
        return (grid(c, i + 1, j, k) - grid(c, i - 1, j, k)) / (2 * grid.dx);
    }

    template<>
    Real d_dy(const Grid<STANDARD> &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j + 1, k) - grid(c, i, j - 1, k)) / (2 * grid.dy);
    }

    template<>
    Real d_dz(const Grid<STANDARD> &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j, k + 1) - grid(c, i, j, k - 1)) / (2 * grid.dz);
    }

    template<>
    Real d2_dx2(const Grid<STANDARD> &grid, Component c, int i, int j, int k) {
        return (grid(c, i - 1, j, k) - 2 * grid(c, i, j, k) + grid(c, i + 1, j, k)) /
               (grid.dx * grid.dx);
    }

    template<>
    Real d2_dy2(const Grid<STANDARD> &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j - 1, k) - 2 * grid(c, i, j, k) + grid(c, i, j + 1, k)) /
               (grid.dy * grid.dy);
    }

    template<>
    Real d2_dz2(const Grid<STANDARD> &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j, k - 1) - 2 * grid(c, i, j, k) + grid(c, i, j, k + 1)) /
               (grid.dz * grid.dz);
    }

    template<>
    Real lap(const Grid<STANDARD> &grid, Component c, int i, int j, int k) {
        return d2_dx2(grid, c, i, j, k) + d2_dy2(grid, c, i, j, k) + d2_dz2(grid, c, i, j, k);
    }

    template<>
    Vector lap(const Grid<STANDARD> &grid, int i, int j, int k) {
        return {
                d2_dx2(grid, U, i, j, k) + d2_dy2(grid, U, i, j, k) + d2_dz2(grid, U, i, j, k),
                d2_dx2(grid, V, i, j, k) + d2_dy2(grid, V, i, j, k) + d2_dz2(grid, V, i, j, k),
                d2_dx2(grid, W, i, j, k) + d2_dy2(grid, W, i, j, k) + d2_dz2(grid, W, i, j, k)
        };
    }

    template<>
    Real
    get_interpolation(const Grid<STANDARD> &grid, Component to, Component from, int i, int j,
                      int k) {
        switch (to) {
            case Component::U:
                if (from == Component::V) {
                    // Interpolate from V grid to U grid
                    return (grid(Component::V, i, j, k) + grid(Component::V, i + 1, j, k) +
                            grid(Component::V, i, j - 1, k) + grid(Component::V, i + 1, j - 1, k)) / 4;
                } else {
                    // Interpolate from W grid to U grid
                    return (grid(Component::W, i, j, k) + grid(Component::W, i + 1, j, k) +
                            grid(Component::W, i, j, k - 1) + grid(Component::W, i + 1, j, k - 1)) / 4;
                }
            case Component::V:
                if (from == Component::U) {
                    // Interpolate from U grid to V grid
                    return (grid(Component::U, i, j, k) + grid(Component::U, i, j + 1, k) +
                            grid(Component::U, i - 1, j, k) + grid(Component::U, i - 1, j + 1, k)) / 4;
                } else {
                    // Interpolate from W grid to V grid
                    return (grid(Component::W, i, j, k) + grid(Component::W, i, j + 1, k) +
                            grid(Component::W, i, j, k - 1) + grid(Component::W, i, j + 1, k - 1)) / 4;
                }
            case Component::W:
                if (from == Component::U) {
                    // Interpolate from U grid to W grid
                    return (grid(Component::U, i, j, k) + grid(Component::U, i, j, k + 1) +
                            grid(Component::U, i - 1, j, k) + grid(Component::U, i - 1, j, k + 1)) / 4;
                } else {
                    // Interpolate from V grid to W grid
                    return (grid(Component::V, i, j, k) + grid(Component::V, i, j, k + 1) +
                            grid(Component::V, i, j - 1, k) + grid(Component::V, i, j - 1, k + 1)) / 4;
                }
            default:
                return 0;
        }
    }

    template<>
    Real conv(const Grid<STANDARD> &grid, Component c, int i, int j, int k) {
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

    template<>
    Vector conv(const Grid<STANDARD> &grid, int i, int j, int k) {
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
}