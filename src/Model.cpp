
#include "Model.hpp"

/**
 * Initialize the space grid with standard addressing
 */
template<>
void Model<STANDARD>::initGrid(const VectorFunction &initial_velocity, const Function &initial_pressure) {
    Real sdx = dx/2;
    Real sdy = dy/2;
    Real sdz = dz/2;

    for (index_t z = 0; z < grid.nz; ++z) {
        for (index_t y = 0; y < grid.ny; ++y) {
            for (index_t x = 0; x < grid.nx; ++x) {
                grid(U, x, y, z) = initial_velocity(x + sdx, y, z)[0];
                grid(V, x, y, z) = initial_velocity(x, y + sdy, z)[1];
                grid(W, x, y, z) = initial_velocity(x, y, z + sdz)[2];
                grid(P, x, y, z) = initial_pressure(x, y, z);
            }
        }
    }
}