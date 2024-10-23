#include <utils.hpp>

using namespace utils;

namespace utils {

    template<>
    Real d_dx(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k) {
        return (model.grid(c, i + 1, j, k) - model.grid(c, i - 1, j, k)) / (2 * model.dx);
    }

    template<>
    Real d_dy(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k) {
        return (model.grid(c, i, j + 1, k) - model.grid(c, i, j - 1, k)) / (2 * model.dy);
    }

    template<>
    Real d_dz(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k) {
        return (model.grid(c, i, j, k + 1) - model.grid(c, i, j, k - 1)) / (2 * model.dz);
    }

    template<>
    Real d2_dx2(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k) {
        return (model.grid(c, i - 1, j, k) - 2 * model.grid(c, i, j, k) + model.grid(c, i + 1, j, k)) /
               (model.dx * model.dx);
    }

    template<>
    Real d2_dy2(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k) {
        return (model.grid(c, i, j - 1, k) - 2 * model.grid(c, i, j, k) + model.grid(c, i, j + 1, k)) /
               (model.dy * model.dy);
    }

    template<>
    Real d2_dz2(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k) {
        return (model.grid(c, i, j, k - 1) - 2 * model.grid(c, i, j, k) + model.grid(c, i, j, k + 1)) /
               (model.dz * model.dz);
    }

    template<>
    Real lap(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k) {
        return d2_dx2(model, c, i, j, k) + d2_dy2(model, c, i, j, k) + d2_dz2(model, c, i, j, k);
    }

    template<Component to, Component from>
    Real get_interpolation(const StaggeredGrid<Addressing_T::STANDARD> &grid, int i, int j, int k) {
        if constexpr (to == U) {
            if (from == V) {
                // Interpolate from V grid to U grid
                return (grid(V, i, j, k) + grid(V, i + 1, j, k) +
                        grid(V, i, j - 1, k) + grid(V, i + 1, j - 1, k)) / 4;
            } else {
                // Interpolate from W grid to U grid
                return (grid(W, i, j, k) + grid(W, i + 1, j, k) +
                        grid(W, i, j, k - 1) + grid(W, i + 1, j, k - 1)) / 4;
            }
        } else if constexpr (to == V) {
            if (from == U) {
                // Interpolate from U grid to V grid
                return (grid(U, i, j, k) + grid(U, i, j + 1, k) +
                        grid(U, i - 1, j, k) + grid(U, i - 1, j + 1, k)) / 4;
            } else {
                // Interpolate from W grid to V grid
                return (grid(W, i, j, k) + grid(W, i, j + 1, k) +
                        grid(W, i, j, k - 1) + grid(W, i, j + 1, k - 1)) / 4;
            }
        } else if constexpr (to == W) {
            if (from == U) {
                // Interpolate from U grid to W grid
                return (grid(U, i, j, k) + grid(U, i, j, k + 1) +
                        grid(U, i - 1, j, k) + grid(U, i - 1, j, k + 1)) / 4;
            } else {
                // Interpolate from V grid to W grid
                return (grid(V, i, j, k) + grid(V, i, j, k + 1) +
                        grid(V, i, j - 1, k) + grid(V, i, j - 1, k + 1)) / 4;
            }
        } else {
            return 0;
        }
    }

    template<Component c>
    Real conv(const Model<Addressing_T::STANDARD> &model, int i, int j, int k) {
        if constexpr (c == U) {
            return model.grid(c, i, j, k) * d_dx(model, c, i, j, k) +
                   get_interpolation<U, V>(model.grid, i, j, k) * d_dy(model, c, i, j, k) +
                   get_interpolation<U, W>(model.grid, i, j, k) * d_dz(model, c, i, j, k);
        } else if constexpr (c == V) {
            return get_interpolation<V, U>(model.grid, i, j, k) * d_dx(model, c, i, j, k) +
                   model.grid(c, i, j, k) * d_dy(model, c, i, j, k) +
                   get_interpolation<V, W>(model.grid, i, j, k) * d_dz(model, c, i, j, k);
        } else if constexpr (c == W) {
            return get_interpolation<W, U>(model.grid, i, j, k) * d_dx(model, c, i, j, k) +
                   get_interpolation<W, V>(model.grid, i, j, k) * d_dy(model, c, i, j, k) +
                   model.grid(c, i, j, k) * d_dz(model, c, i, j, k);
        } else {
            return 0;
        }
    }
}