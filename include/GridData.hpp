#ifndef AEROHPC_A_GRIDDATA_H
#define AEROHPC_A_GRIDDATA_H

#include <vector>
#include <cstdio>
#include <array>
#include "Traits.hpp"
#include "GridStructure.hpp"

class GridData {

    /**
     * Nodes data
     */
    Real *u;
    Real *v;
    Real *w;
    Real *p;

    /**
     * Calculate the indexing used into velocity and pressure arrays
     */
    index_t indexing(index_t x, index_t y, index_t z) const;

public:

    GridData() = delete;

    ~GridData() {
        delete[] u;
        delete[] v;
        delete[] w;
        delete[] p;
    }

    /**
     * Construct a staggered grid with ghost nodes,
     * @param nodes number of non-ghost nodes
     * @param spacing dx, dy, dz info
     * @param ghosts number of ghost nodes
     */
    explicit GridData(const GridStructure& structure) : structure(structure){
        auto dim = structure.grid_nodes[0] * structure.grid_nodes[1] * structure.grid_nodes[2];
        u = new Real[dim];
        v = new Real[dim];
        w = new Real[dim];
        p = new Real[dim];
    }

    /**
     * Reference Grid Structure
     */
     const GridStructure &structure;

    /**
     * Operator that accesses the memory using a 3D view of the object
     */
#define get_component(name, comp) \
inline Real &name(index_t i, index_t j, index_t k){ \
    return comp[indexing(i,j,k)]; \
} \
inline Real &name(index_t i, index_t j, index_t k) const { \
    return comp[indexing(i,j,k)]; \
}

    get_component(U, u)

    get_component(V, v)

    get_component(W, w)

    get_component(P, p)

    /**
     * Initialize the grid given the initial velocity and pressure function
     */
    void initData(const VectorFunction &initial_velocity, const Function &initial_pressure);

    /**
     * Swap grid data with the given one
     */
    void swap(GridData &other) noexcept {
        std::swap(u, other.u);
        std::swap(v, other.v);
        std::swap(w, other.w);
        std::swap(p, other.p);
    }
};

#endif //AEROHPC_A_GRIDDATA_H
