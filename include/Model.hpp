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
    std::vector<BoundaryCondition<A> *> _bcs;

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
          const VectorFunction &initial_velocity, const Function &initial_pressure)
            : grid(grid),
              reynolds(reynolds),
              spacing(spacing),
              staggered_spacing{spacing[0] / 2,
                                spacing[1] / 2,
                                spacing[2] / 2} {
        initGrid(initial_velocity, initial_pressure);
    }

    Model(const Vector &spacing, StaggeredGrid<A> &grid, Real reynolds)
            : grid(grid),
              reynolds(reynolds),
              spacing(spacing),
              staggered_spacing{spacing[0] / 2,
                                spacing[1] / 2,
                                spacing[2] / 2} {}
    /**
     * Add a boundary condition to the list of the model
     */
    void addBC(BoundaryCondition<A> &bc) { _bcs.push_back(&bc); }

    /**
     * Apply all the boundary conditions of the model,
     * Be aware that the insertion order of the conditions can affect the result of the application
     */
    void applyBCs(Real time) { for (BoundaryCondition<A> *bc: _bcs) bc->apply(grid, time); }

    /**
      * Model space grid
      */
    StaggeredGrid<A> &grid;

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
     * Node spacing
     */
    const Vector staggered_spacing;

    /**
     * Alias for x-axes spacing
     */
    const Real &sdx = staggered_spacing[0];

    /**
     * Alias for y-axes spacing
     */
    const Real &sdy = staggered_spacing[1];

    /**
     * Alias for z-axes spacing
     */
    const Real &sdz = staggered_spacing[2];

    /**
     * Reynolds number
     */
    const Real reynolds;
};

#endif //AEROHPC_A_MODEL_H
