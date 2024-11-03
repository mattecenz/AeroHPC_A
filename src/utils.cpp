#include <utils.hpp>

using namespace utils;

namespace utils {
    Real d_dx(const Grid &grid, Component c, int i, int j, int k) {
        return (grid(c, i + 1, j, k) - grid(c, i - 1, j, k)) / (2 * grid.dx);
    }

    Real d_dy(const Grid &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j + 1, k) - grid(c, i, j - 1, k)) / (2 * grid.dy);
    }

    Real d_dz(const Grid &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j, k + 1) - grid(c, i, j, k - 1)) / (2 * grid.dz);
    }

    Real d2_dx2(const Grid &grid, Component c, int i, int j, int k) {
        return (grid(c, i - 1, j, k) - 2 * grid(c, i, j, k) + grid(c, i + 1, j, k)) /
               (grid.dx * grid.dx);
    }

    Real d2_dy2(const Grid &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j - 1, k) - 2 * grid(c, i, j, k) + grid(c, i, j + 1, k)) /
               (grid.dy * grid.dy);
    }

    Real d2_dz2(const Grid &grid, Component c, int i, int j, int k) {
        return (grid(c, i, j, k - 1) - 2 * grid(c, i, j, k) + grid(c, i, j, k + 1)) /
               (grid.dz * grid.dz);
    }

    // template<>
    // Real lap(const Grid &grid, Component c, int i, int j, int k) {
    //     return d2_dx2(grid, c, i, j, k) + d2_dy2(grid, c, i, j, k) + d2_dz2(grid, c, i, j, k);
    // }

    Vector lap(const Grid &grid, int i, int j, int k) {
        return {
            d2_dx2(grid, U, i, j, k) + d2_dy2(grid, U, i, j, k) + d2_dz2(grid, U, i, j, k),
            d2_dx2(grid, V, i, j, k) + d2_dy2(grid, V, i, j, k) + d2_dz2(grid, V, i, j, k),
            d2_dx2(grid, W, i, j, k) + d2_dy2(grid, W, i, j, k) + d2_dz2(grid, W, i, j, k)
        };
    }

    template <>
    Real lap<U>(const Grid &grid, const int i, const int j, const int k) {
        return d2_dx2(grid, U, i, j, k) + d2_dy2(grid, U, i, j, k) + d2_dz2(grid, U, i, j, k);
    }



    Real lap_u(const Grid &grid, const int i, const int j, const int k) {
        return d2_dx2(grid, U, i, j, k) + d2_dy2(grid, U, i, j, k) + d2_dz2(grid, U, i, j, k);
    }


    template<>
    Real lap<V>(const Grid &grid, const int i, const int j, const int k) {
        return d2_dx2(grid, V, i, j, k) + d2_dy2(grid, V, i, j, k) + d2_dz2(grid, V, i, j, k);
    }

    template<>
    Real lap<W>(const Grid &grid, const int i, const int j, const int k) {
        return d2_dx2(grid, W, i, j, k) + d2_dy2(grid, W, i, j, k) + d2_dz2(grid, W, i, j, k);
    }


    template<Component to, Component from>
    Real get_interpolation(const Grid &grid, const int i, const int j, const int k) {
        if constexpr (to == U) {
            if constexpr (from == V) {
                return (grid(Component::V, i, j, k) + grid(Component::V, i + 1, j, k) +
                        grid(Component::V, i, j - 1, k) + grid(Component::V, i + 1, j - 1, k)) / 4;
            } else {
                return (grid(Component::W, i, j, k) + grid(Component::W, i + 1, j, k) +
                        grid(Component::W, i, j, k - 1) + grid(Component::W, i + 1, j, k - 1)) / 4;
            }
        } else if constexpr (to == V) {
            if constexpr (from == U) {
                return (grid(Component::U, i, j, k) + grid(Component::U, i, j + 1, k) +
                        grid(Component::U, i - 1, j, k) + grid(Component::U, i - 1, j + 1, k)) / 4;
            } else {
                return (grid(Component::W, i, j, k) + grid(Component::W, i, j + 1, k) +
                        grid(Component::W, i, j, k - 1) + grid(Component::W, i, j + 1, k - 1)) / 4;
            }
        } else {
            if constexpr (from == U) {
                return (grid(Component::U, i, j, k) + grid(Component::U, i, j, k + 1) +
                        grid(Component::U, i - 1, j, k) + grid(Component::U, i - 1, j, k + 1)) / 4;
            } else {
                return (grid(Component::V, i, j, k) + grid(Component::V, i, j, k + 1) +
                        grid(Component::V, i, j - 1, k) + grid(Component::V, i, j - 1, k + 1)) / 4;
            }
        }
    }

    template<>
    Real conv<U>(const Grid &grid, const int i, const int j, const int k) {
        return grid(U, i, j, k) * d_dx(grid, U, i, j, k) +
               get_interpolation<U, V>(grid, i, j, k) * d_dy(grid, U, i, j, k) +
               get_interpolation<U, W>(grid, i, j, k) * d_dz(grid, U, i, j, k);
    }

    template<>
    Real conv<V>(const Grid &grid, const int i, const int j, const int k) {
        return get_interpolation<V, U>(grid, i, j, k) * d_dx(grid, V, i, j, k) +
               grid(V, i, j, k) * d_dy(grid, V, i, j, k) +
               get_interpolation<V, W>(grid, i, j, k) * d_dz(grid, V, i, j, k);
    }

    template<>
    Real conv<W>(const Grid &grid, const int i, const int j, const int k) {
        return get_interpolation<W, U>(grid, i, j, k) * d_dx(grid, W, i, j, k) +
               get_interpolation<W, V>(grid, i, j, k) * d_dy(grid, W, i, j, k) +
               grid(W, i, j, k) * d_dz(grid, W, i, j, k);
    }
}
