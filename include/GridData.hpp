#ifndef AEROHPC_A_GRIDDATA_H
#define AEROHPC_A_GRIDDATA_H

#include "Traits.hpp"
#include "GridStructure.hpp"

class GridData {

    /**
     * Calculate the indexing used into velocity and pressure arrays
     */
    index_t indexing(index_t x, index_t y, index_t z) const;


    /**
     * Dimension of one component of data considering the grid nodes with ghosts
     */
    const index_t grid_dim;

public:

    /**
     * Dimension of one component of data considering the grid nodes without ghosts
     */
    const index_t node_dim;

    /**
     * Data (concatenation of flattened array of all velocity components + optional pressure)
     */
    Real *data;

    GridData() = delete;

    /**
     * Construct a staggered grid with ghost nodes,
     * @param nodes number of non-ghost nodes
     * @param spacing dx, dy, dz info
     * @param ghosts number of ghost nodes
     */
    GridData(const GridStructure &structure, bool allocPressure) : structure(structure),
                                                                   node_dim(structure.nx *
                                                                            structure.ny *
                                                                            structure.nz),
                                                                   grid_dim(structure.gx *
                                                                            structure.gy *
                                                                            structure.gz) {
        data = new Real[grid_dim * (allocPressure ? 4 : 3)];
    }

    explicit GridData(const GridStructure &structure) : GridData(structure, true) {}

    /**
     * Grid Structure
     */
    const GridStructure structure;

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

    get_component(U, data)

    get_component(V, (&data[grid_dim]))

    get_component(W, (&data[grid_dim * 2]))

    get_component(P, (&data[grid_dim * 3]))

    /**
     * Initialize the grid given the initial velocity and pressure function
     */
    void initData(const VectorFunction &initial_velocity, const Function &initial_pressure);

    /**
     * Swap grid data with the given one
     */
    void swap(GridData &other) noexcept {
        std::swap(data, other.data);
    }
};

#endif //AEROHPC_A_GRIDDATA_H
