#ifndef AEROHPC_A_STAGGEREDGRID_H
#define AEROHPC_A_STAGGEREDGRID_H


#include <vector>
#include <cstdio>
#include <array>

#include <Traits.hpp>

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

/**
 *
 * Class for handling a 3D Matrix addressing.
 *
 * Used as layer of abstraction in order to separate the addressing complexity
 * from the maths development.
 *
 * For the moment the ordering is: first x-y plane as if it were a row-major matrix then move on the z axis
 *
 */
template<Addressing_T Addressing>
class StaggeredGrid {

public:

    // Constructors

    StaggeredGrid() = delete;

    /**
     * Constructor of a nx * ny * nz tensor with initial value equal to zero
     */
    explicit StaggeredGrid(const std::array<size_t, 3> &nodes) :
            nodes(nodes) {
        auto dim = nodes[0] * nodes[1] * nodes[2] * Component::N_COMPONENTS;
        _entries.resize(dim);
        std::fill(_entries.begin(), _entries.end(), 0);
    }

#pragma inline
    /**
     * Operator that accesses the memory using a 3D view of the object
     */
    Real &operator()(Component c, size_t i, size_t j, size_t k);

#pragma inline
    /**
     * Operator that accesses 3D view for read-only operations
     */
    const Real &operator()(Component c, size_t i, size_t j, size_t k) const;

    /**
     * Number of nodes
     */
    const std::array<size_t, 3> nodes;

    /**
     * Alias for x-axes nodes number
     */
    const size_t &nx = nodes[0];
    /**
     * Alias for y-axes nodes number
     */
    const size_t &ny = nodes[1];
    /**
     * Alias for z-axes nodes number
     */
    const size_t &nz = nodes[2];

private:
    /**
     * Nodes data
     */
    std::vector<Real> _entries;

};

#endif //AEROHPC_A_STAGGEREDGRID_H
