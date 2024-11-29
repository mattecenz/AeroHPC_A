#ifndef AEROHPC_A_GRIDSTRUCTURE_H
#define AEROHPC_A_GRIDSTRUCTURE_H

#include "Traits.hpp"

class GridStructure {

public:
    GridStructure(const Idx3 &nodes, const std::array<Real, 3> &spacing, const Idx3 &displacement, const index_t ghosts)
            : nodes(nodes),
              gp(ghosts),
              grid_nodes({nodes[0] + (2 * ghosts),
                           nodes[1] + (2 * ghosts),
                           nodes[2] + (2 * ghosts)}),
              spacing(spacing),
              displacement(displacement),
              staggered_spacing{spacing[0] / real(2),
                                spacing[1] / real(2),
                                spacing[2] / real(2)} {}

    /**
    * Padding induced by ghost point
    */
    const index_t gp;

    /**
     * Number of total nodes (with ghosts)
     */
    const Idx3 grid_nodes;
    /**
     * Alias for x-axes total nodes number
     */
    const index_t &gx = grid_nodes[0];
    /**
     * Alias for y-axes total nodes number
     */
    const index_t &gy = grid_nodes[1];
    /**
     * Alias for z-axes total nodes number
     */
    const index_t &gz = grid_nodes[2];

    /**
     * Number of nodes (without ghost points)
     */
    const Idx3 nodes;
    /**
     * Alias for x-axes nodes number
     */
    const index_t &nx = nodes[0];
    /**
     * Alias for y-axes nodes number
     */
    const index_t &ny = nodes[1];
    /**
     * Alias for z-axes nodes number
     */
    const index_t &nz = nodes[2];

    /**
     * Global displacement of the grid
     */
    const Idx3 displacement;
    /**
     * Alias for x-axes nodes number
     */
    const index_t &px = displacement[0];
    /**
     * Alias for y-axes nodes number
     */
    const index_t &py = displacement[1];
    /**
     * Alias for z-axes nodes number
     */
    const index_t &pz = displacement[2];

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
};

#endif //AEROHPC_A_GRIDSTRUCTURE_H
