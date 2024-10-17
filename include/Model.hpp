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
    Model(const std::array<Real, 3>& phy_dim, const std::array<size_t, 3>& nodes, Real Reynolds,
          const Function &initial_velocity, const Function &initial_pressure,
          Function boundary_inlet) : _boundary_inlet(std::move(boundary_inlet)),
                                     _Reynolds(Reynolds),
                                     _grid(nodes),
                                     _spacing{phy_dim[0] / nodes[0],
                                              phy_dim[1] / nodes[1],
                                              phy_dim[2] / nodes[2]} {
        initGrid(initial_velocity, initial_pressure);
    }

    /**
     * Returns the grid of the model
     */
    StaggeredGrid<A> &grid() { return _grid; }

    /**
     * Returns the node spacing of the model
     */
    const std::array<Real, 3> &spacing() { return _spacing; }

private:
    /**
     * Node spacing
     */
    const std::array<Real, 3> _spacing;

    /**
     * Reynolds number
     */
    const Real _Reynolds;
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
