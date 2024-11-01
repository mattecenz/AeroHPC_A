#ifndef AEROHPC_A_GHOSTEDSG_H
#define AEROHPC_A_GHOSTEDSG_H

#include <vector>
#include <cstdio>
#include <array>

#include "Traits.hpp"

/**
 * Components of the problem
 */
enum Component {
    U = 0, // x-axis velocity
    V = 1, // y-axis velocity
    W = 2, // z-axis velocity
    P = 3, // pressure
    N_COMPONENTS = 4
};


template<Addressing_T A>
class Grid {

private:
    /**
     * Padding induced by ghost point
     */
    index_t _gp;

    /**
     * Nodes data
     */
    std::vector<Real> _entries;

    /**
     * Number of nodes without ghost points
     */
    std::array<index_t, 3> _nodes;

    /**
     * Number of nodes with ghost points
     */
    std::array<index_t, 3> _grid_nodes;


public:

    Grid() = delete;

    /**
     * Construct a staggered grid with ghost nodes,
     * @param nodes number of non-ghost nodes
     * @param ghosts_p_axis_p_direction number of ghost nodes,
     *                                  3 couples (one per axis),
     *                                  each couple define the number of ghost nodes (before, after) non-ghost nodes
     */
    Grid(const std::array<index_t, 3> &nodes, const std::array<Real,3> spacing, const index_t ghosts) :
            _nodes(nodes),
            _gp(ghosts),
            _grid_nodes({nodes[0] + (2 * ghosts),
                         nodes[1] + (2 * ghosts),
                         nodes[2] + (2 * ghosts)}),
            spacing(spacing),
            staggered_spacing{spacing[0] / real(2),
                              spacing[1] / real(2),
                              spacing[2] / real(2)} {
        auto dim = _grid_nodes[0] * _grid_nodes[1] * _grid_nodes[2] * N_COMPONENTS;
        _entries.resize(dim);
    };

    /**
     * Operator that accesses the memory using a 3D view of the object
     */
    inline Real &operator()(Component c, index_t i, index_t j, index_t k){
        return _entries[c + ((i + _gp) + ((j + _gp) + (k + _gp) * _grid_nodes[1]) * _grid_nodes[0]) * N_COMPONENTS];
    }


    /**
     * Operator that accesses 3D view for read-only operations
     */
    inline const Real &operator()(Component c, index_t i, index_t j, index_t k) const {
        return _entries[c + ((i + _gp) + ((j + _gp) + (k + _gp) * _grid_nodes[1]) * _grid_nodes[0]) * N_COMPONENTS];
    }

    /**
     * Initialize the grid given the initial velocity and pressure function
     */
    void initGrid(const VectorFunction &initial_velocity, const Function &initial_pressure);

    /**
     * Number of nodes
     */
    const std::array<index_t, 3> &nodes = _nodes;

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
     * Alias for ghost point number
     */
    const index_t &gp = _gp;

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

#endif //AEROHPC_A_GHOSTEDSG_H
