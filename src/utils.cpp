#include <utils.hpp>
using namespace utils;

namespace utils
{

    template <>
    Real d_dx(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k)
    {
        return (model.grid(c, i + 1, j, k) - model.grid(c, i - 1, j, k)) / (2*model.dx);
    }

    template <>
    Real d_dy(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k)
    {
        return (model.grid(c, i, j + 1, k) - model.grid(c, i, j - 1, k)) / (2*model.dy);
    }

    template <>
    Real d_dz(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k)
    {
        return (model.grid(c, i, j, k + 1) - model.grid(c, i, j, k - 1)) / (2*model.dz);
    }

    template <>
    Real d2_dx2(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k)
    {
        return (model.grid(c, i - 1, j, k) - 2 * model.grid(c, i, j, k) + model.grid(c, i + 1, j, k)) / (model.dx * model.dx);
    }

    template <>
    Real d2_dy2(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k)
    {
        return (model.grid(c, i, j - 1, k) - 2 * model.grid(c, i, j, k) + model.grid(c, i, j + 1, k)) / (model.dy * model.dy);
    }

    template <>
    Real d2_dz2(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k)
    {
        return (model.grid(c, i, j, k - 1) - 2 * model.grid(c, i, j, k) + model.grid(c, i, j, k + 1)) / (model.dz * model.dz);
    }

    template <>
    Real lap(const Model<Addressing_T::STANDARD> &model, Component c, int i,int j,int k){
        return d2_dx2(model,c,i,j,k)+d2_dy2(model,c,i,j,k)+d2_dz2(model,c,i,j,k);
    }

    template <>
    Real get_interpolation(const StaggeredGrid<Addressing_T::STANDARD> &grid, Component to, Component from, int i, int j, int k)
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

    template <>
    Real conv(const Model<Addressing_T::STANDARD> &model, Component c, int i, int j, int k)
    {
        switch (c)
        {
        case Component::U:
            return model.grid(c, i, j, k) * d_dx(model, c, i, j, k) +
                get_interpolation(model.grid, Component::U, Component::V, i, j, k) * d_dy(model, c, i, j, k) +
                get_interpolation(model.grid, Component::U, Component::W, i, j, k) * d_dz(model, c, i, j, k);
            break;
        case Component::V:
            return get_interpolation(model.grid, Component::V, Component::U, i, j, k) * d_dx(model, c, i, j, k) +
                model.grid(c, i, j, k) * d_dy(model, c, i, j, k) +
                get_interpolation(model.grid, Component::V, Component::W, i, j, k) * d_dz(model, c, i, j, k);
            break;
        case Component::W:
            return get_interpolation(model.grid, Component::W, Component::U, i, j, k) * d_dx(model, c, i, j, k) +
                get_interpolation(model.grid, Component::W, Component::V, i, j, k) * d_dy(model, c, i, j, k) +
                model.grid(c, i, j, k) * d_dz(model, c, i, j, k);
            break;
        }
    }
}