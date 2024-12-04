#ifndef AEROHPC_A_GRIDDATA_H
#define AEROHPC_A_GRIDDATA_H

#include "Traits.hpp"
#include "GridStructure.hpp"

class GridData {

    /**
     * Calculate the indexing used into velocity and pressure arrays
     */
    index_t indexing(index_t x, index_t y, index_t z) const;

public:

    /**
     * Dimension of one component of data
     */
    const index_t dim;

    /**
     * Velocity Data (concatenation of flattened array of all velocity components)
     */
    Real *velocity_data;

    /**
     * Pressure Data (flattened array of pressure)
     */
    Real *pressure_data;


    GridData() = delete;

    ~GridData() {
        delete[] velocity_data;
    }

    /**
     * Construct a staggered grid with ghost nodes,
     * @param nodes number of non-ghost nodes
     * @param spacing dx, dy, dz info
     * @param ghosts number of ghost nodes
     */
    GridData(const GridStructure &structure, bool allocPressure) : structure(structure),
                                                                            dim(structure.grid_nodes[0] * structure.grid_nodes[1] *
                                                                                structure.grid_nodes[2]) {
        velocity_data = new Real[dim * 3];
        pressure_data = nullptr;
        if (allocPressure)
            pressure_data = new Real[dim];
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

    get_component(U, velocity_data)

    get_component(V, (&velocity_data[dim]))

    get_component(W, (&velocity_data[dim + dim]))

    get_component(P, pressure_data)

    /**
     * Initialize the grid given the initial velocity and pressure function
     */
    void initData(const VectorFunction &initial_velocity, const Function &initial_pressure);

    /**
     * Swap grid data with the given one
     */
    void swap(GridData &other) noexcept {
        std::swap(velocity_data, other.velocity_data);
        std::swap(pressure_data, other.pressure_data);
    }
};

#endif //AEROHPC_A_GRIDDATA_H
