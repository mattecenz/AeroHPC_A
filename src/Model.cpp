
#include "Model.hpp"

/**
 * Initialize the space grid with standard addressing
 */
template<>
void Model<STANDARD>::initGrid(const VectorFunction &initial_velocity, const Function &initial_pressure) {
    for (index_t z = 0; z < grid.nz; ++z) {
        for (index_t y = 0; y < grid.ny; ++y) {
            for (index_t x = 0; x < grid.nx; ++x) {
                Real px = real(x) * dx;
                Real py = real(y) * dy;
                Real pz = real(z) * dz;

                grid(U, x, y, z) = initial_velocity(px + dx, py + sdy, pz + sdz)[0];
                grid(V, x, y, z) = initial_velocity(px + sdx, py + dy, pz + sdz)[1];
                grid(W, x, y, z) = initial_velocity(px + sdx, py + sdy, pz + dz)[2];
                grid(P, x, y, z) = initial_pressure(px + sdx, py + sdy, pz + sdz);
            }
        }
    }
}