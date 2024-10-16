#ifndef AEROHPC_A_STAGGEREDGRID_H
#define AEROHPC_A_STAGGEREDGRID_H


#include <vector>
#include <cstdio>

#include <Traits.hpp>

/**
 * Components of the problem
 */
enum Component{
    U = 0,
    V = 1,
    W = 2,
    P = 3,
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
template<typename Type, Addressing_T Addressing>
class StaggeredGrid {

public:

    // Constructors

    /**
     * Constructor of a nx * ny * nz tensor with initial value equal to zero
     */
    StaggeredGrid(const size_t nx, const size_t ny, const size_t nz) :
            _nx(nx),
            _ny(ny),
            _nz(nz) {
        auto dim = nx * ny * nz * Component::N_COMPONENTS;
        _entries.resize(dim);
        std::fill(_entries.begin(), _entries.end(), 0);

    }

    /**
     * Operator that accesses the memory using a 3D view of the object
     */
    Type &operator()(Component c, size_t i, size_t j, size_t k);

    /**
     * Operator that accesses 3D view for read-only operations
     */
    const Type &operator()(Component c, size_t i, size_t j, size_t k) const;

private:
    /**
     * Nodes data
     */
    std::vector<Real> _entries;

    /**
     * Number of nodes on x axis
     */
    size_t _nx;
    /**
     * Number of nodes on y axis
     */
    size_t _ny;
    /**
     * Number of nodes on z axis
     */
    size_t _nz;
};


#endif //AEROHPC_A_STAGGEREDGRID_H
