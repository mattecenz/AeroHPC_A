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

    // template<>
    // Real lap(const Grid<STANDARD> &grid, Component c, int i, int j, int k) {
    //     return d2_dx2(grid, c, i, j, k) + d2_dy2(grid, c, i, j, k) + d2_dz2(grid, c, i, j, k);
    // }

    template<>
    Vector lap(const Grid<STANDARD> &grid, int i, int j, int k) {
        return {
            d2_dx2(grid, U, i, j, k) + d2_dy2(grid, U, i, j, k) + d2_dz2(grid, U, i, j, k),
            d2_dx2(grid, V, i, j, k) + d2_dy2(grid, V, i, j, k) + d2_dz2(grid, V, i, j, k),
            d2_dx2(grid, W, i, j, k) + d2_dy2(grid, W, i, j, k) + d2_dz2(grid, W, i, j, k)
        };
    }

    template<>
    Real lap<U>(const Grid<STANDARD> &grid, const int i, const int j, const int k) {
        return d2_dx2(grid, U, i, j, k) + d2_dy2(grid, U, i, j, k) + d2_dz2(grid, U, i, j, k);
    }



    Real lap_u(const Grid<STANDARD> &grid, const int i, const int j, const int k) {
        return d2_dx2(grid, U, i, j, k) + d2_dy2(grid, U, i, j, k) + d2_dz2(grid, U, i, j, k);
    }


    template<>
    Real lap<V>(const Grid<STANDARD> &grid, const int i, const int j, const int k) {
        return d2_dx2(grid, V, i, j, k) + d2_dy2(grid, V, i, j, k) + d2_dz2(grid, V, i, j, k);
    }

    template<>
    Real lap<W>(const Grid<STANDARD> &grid, const int i, const int j, const int k) {
        return d2_dx2(grid, W, i, j, k) + d2_dy2(grid, W, i, j, k) + d2_dz2(grid, W, i, j, k);
    }


    template<Component to, Component from>
    Real get_interpolation(const Grid<STANDARD> &grid, const int i, const int j, const int k) {
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

    // template<>
    // Real
    // get_interpolation(const Grid<STANDARD> &grid, Component to, Component from, int i, int j,
    //                   int k) {
    //     switch (to) {
    //         case Component::U:
    //             if (from == Component::V) {
    //                 // Interpolate from V grid to U grid
    //                 return (grid(Component::V, i, j, k) + grid(Component::V, i + 1, j, k) +
    //                         grid(Component::V, i, j - 1, k) + grid(Component::V, i + 1, j - 1, k)) / 4;
    //             } else {
    //                 // Interpolate from W grid to U grid
    //                 return (grid(Component::W, i, j, k) + grid(Component::W, i + 1, j, k) +
    //                         grid(Component::W, i, j, k - 1) + grid(Component::W, i + 1, j, k - 1)) / 4;
    //             }
    //         case Component::V:
    //             if (from == Component::U) {
    //                 // Interpolate from U grid to V grid
    //                 return (grid(Component::U, i, j, k) + grid(Component::U, i, j + 1, k) +
    //                         grid(Component::U, i - 1, j, k) + grid(Component::U, i - 1, j + 1, k)) / 4;
    //             } else {
    //                 // Interpolate from W grid to V grid
    //                 return (grid(Component::W, i, j, k) + grid(Component::W, i, j + 1, k) +
    //                         grid(Component::W, i, j, k - 1) + grid(Component::W, i, j + 1, k - 1)) / 4;
    //             }
    //         case Component::W:
    //             if (from == Component::U) {
    //                 // Interpolate from U grid to W grid
    //                 return (grid(Component::U, i, j, k) + grid(Component::U, i, j, k + 1) +
    //                         grid(Component::U, i - 1, j, k) + grid(Component::U, i - 1, j, k + 1)) / 4;
    //             } else {
    //                 // Interpolate from V grid to W grid
    //                 return (grid(Component::V, i, j, k) + grid(Component::V, i, j, k + 1) +
    //                         grid(Component::V, i, j - 1, k) + grid(Component::V, i, j - 1, k + 1)) / 4;
    //             }
    //         default:
    //             return 0;
    //     }
    // }

    // template<>
    // Real conv(const Grid<STANDARD> &grid, Component c, int i, int j, int k) {
    //     switch (c) {
    //         case Component::U:
    //             return grid(c, i, j, k) * d_dx(grid, c, i, j, k) +
    //                    get_interpolation<U, V>(grid, i, j, k) * d_dy(grid, c, i, j, k) +
    //                    get_interpolation<U, W>(grid, i, j, k) * d_dz(grid, c, i, j, k);
    //         case Component::V:
    //             return get_interpolation<V, U>(grid, i, j, k) * d_dx(grid, c, i, j, k) +
    //                    grid(c, i, j, k) * d_dy(grid, c, i, j, k) +
    //                    get_interpolation<V, W>(grid, i, j, k) * d_dz(grid, c, i, j, k);
    //         case Component::W:
    //             return get_interpolation<W, U>(grid, i, j, k) * d_dx(grid, c, i, j, k) +
    //                    get_interpolation<W, V>(grid, i, j, k) * d_dy(grid, c, i, j, k) +
    //                    grid(c, i, j, k) * d_dz(grid, c, i, j, k);
    //         default:
    //             return 0;
    //     }
    // }

    template<>
    Vector conv(const Grid<STANDARD> &grid, int i, int j, int k) {
        return {
            grid(U, i, j, k) * d_dx(grid, U, i, j, k) +
            get_interpolation<U, V>(grid, i, j, k) * d_dy(grid, U, i, j, k) +
            get_interpolation<U, W>(grid, i, j, k) * d_dz(grid, U, i, j, k),
            get_interpolation<V, U>(grid, i, j, k) * d_dx(grid, V, i, j, k) +
            grid(V, i, j, k) * d_dy(grid, V, i, j, k) +
            get_interpolation<V, W>(grid, i, j, k) * d_dz(grid, V, i, j, k),
            get_interpolation<W, U>(grid, i, j, k) * d_dx(grid, W, i, j, k) +
            get_interpolation<W, V>(grid, i, j, k) * d_dy(grid, W, i, j, k) +
            grid(W, i, j, k) * d_dz(grid, W, i, j, k)
        };
    }


    Real conv_u(const Grid<STANDARD> &grid, const int i, const int j, const int k) {
        return grid(U, i, j, k) * d_dx(grid, U, i, j, k) +
               get_interpolation<U, V>(grid, i, j, k) * d_dy(grid, U, i, j, k) +
               get_interpolation<U, W>(grid, i, j, k) * d_dz(grid, U, i, j, k);
    }

    template<>
    Real conv<V>(const Grid<STANDARD> &grid, const int i, const int j, const int k) {
        return get_interpolation<V, U>(grid, i, j, k) * d_dx(grid, V, i, j, k) +
               grid(V, i, j, k) * d_dy(grid, V, i, j, k) +
               get_interpolation<V, W>(grid, i, j, k) * d_dz(grid, V, i, j, k);
    }

    template<>
    Real conv<W>(const Grid<STANDARD> &grid, const int i, const int j, const int k) {
        return get_interpolation<W, U>(grid, i, j, k) * d_dx(grid, W, i, j, k) +
               get_interpolation<W, V>(grid, i, j, k) * d_dy(grid, W, i, j, k) +
               grid(W, i, j, k) * d_dz(grid, W, i, j, k);
    }

    // template<Component C>
    // Real conv(const Grid<STANDARD> &grid, const int i, const int j, const int k) {
    //     if constexpr (C == U) {
    //         return grid(U, i, j, k) * d_dx(grid, U, i, j, k) +
    //                get_interpolation<U, V>(grid, i, j, k) * d_dy(grid, U, i, j, k) +
    //                get_interpolation<U, W>(grid, i, j, k) * d_dz(grid, U, i, j, k);
    //     } else if constexpr (C == V) {
    //         return get_interpolation<V, U>(grid, i, j, k) * d_dx(grid, V, i, j, k) +
    //                grid(V, i, j, k) * d_dy(grid, V, i, j, k) +
    //                get_interpolation<V, W>(grid, i, j, k) * d_dz(grid, V, i, j, k);
    //     } else {
    //         return get_interpolation<W, U>(grid, i, j, k) * d_dx(grid, W, i, j, k) +
    //                get_interpolation<W, V>(grid, i, j, k) * d_dy(grid, W, i, j, k) +
    //                grid(W, i, j, k) * d_dz(grid, W, i, j, k);
    //     }
    // }
}
