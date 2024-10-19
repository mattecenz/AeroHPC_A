#ifndef AEROHPC_A_MODEL_H
#define AEROHPC_A_MODEL_H

#include "Traits.hpp"
#include "StaggeredGrid.hpp"
#include <array>
#include <functional>
#include <utility>

/**
 * Typedef shortening lambda definition
 */
typedef std::function<Real(Real x, Real y, Real z)> Function;

/**
 * Model of the problem
 */
template<Addressing_T A>
class Model {
public:
    /**
     * Construct a model given:
     *  - physical dimension
     *  - number of nodes for each axes
     *  - Reynolds number
     *  - Initial velocity function
     *  - Initial pressure function
     *  - Inlet velocity function
     */
    Model(const std::array<Real, 3> &spacing, const std::array<size_t, 3> &nodes, Real reynolds,
          const Function &initial_velocity, const Function &initial_pressure,
          Function boundary_inlet) : _boundary_inlet(std::move(boundary_inlet)),
                                     reynolds(reynolds),
                                     _grid(nodes),
                                     spacing(spacing) {
        initGrid(initial_velocity, initial_pressure);
    }

    /**
     * Returns the grid of the model
     */
    StaggeredGrid<A> &grid = _grid;

    /**
     * Node spacing
     */
    const std::array<Real, 3> spacing;

    /**
     * Alias for x-axes spacing
     */
    const Real &dx = spacing[0];

    /**
     * Alias for y-axes spacing
     */
    const Real &dy = spacing[1];

    /**
     * Alias for z-axes spacing
     */
    const Real &dz = spacing[2];

    /**
     * Reynolds number
     */
    const Real reynolds;

private:

    /**
     * Inlet velocity boundary condition
     */
    const Function _boundary_inlet;
    /**
     * Model space grid
     */
    StaggeredGrid<A> _grid;

    /**
     * Initialize the space grid
     */
    void initGrid(const Function &initial_velocity, const Function &initial_pressure);
};

#endif //AEROHPC_A_MODEL_H
