#ifndef AEROHPC_A_MODEL_H
#define AEROHPC_A_MODEL_H

#include "Traits.hpp"
#include "StaggeredGrid.hpp"
#include "BoundaryCondition.hpp"
#include <array>
#include <utility>


/**
 * Model of the problem
 */
template<Addressing_T A>
class Model {
private:

    /**
     * Inlet velocity boundary condition
     */
    std::vector<BoundaryCondition<A>> _bcs;

    /**
     * Model space grid
     */
    StaggeredGrid<A> *_grid;

    /**
     * Initialize the space grid
     */
    void initGrid(const VectorFunction &initial_velocity, const Function &initial_pressure);


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
    Model(const Vector &spacing, StaggeredGrid<A> &grid, Real reynolds,
          const VectorFunction &initial_velocity, const Function &initial_pressure) : _grid(&grid),
                                                                                      reynolds(reynolds),
                                                                                      spacing(spacing) {
        initGrid(initial_velocity, initial_pressure);
    }

    // Copy constructor, used in the RK method
    // Leaves the grid empty since it will be overwritten anyway
    Model(Model& m) : _grid(m.grid.nodes), spacing(m.spacing), reynolds(m.reynolds){}

    /**
     * Add a boundary condition to the list of the model
     */
    void addBC(const BoundaryCondition<A> &bc) { _bcs.push_back(bc); }

    /**
     * Apply all the boundary conditions of the model,
     * Be aware that the insertion order of the conditions can affect the result of the application
     */
    void applyBCs() { for (BoundaryCondition<A> bc: _bcs) bc.apply(_grid); }

    /**
     * Returns the grid of the model
     */
    StaggeredGrid<A> grid = *_grid;

    /**
     * Node spacing
     */
    const Vector spacing;

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
};

#endif //AEROHPC_A_MODEL_H
