#include "GridData.hpp"

index_t GridData::indexing(const index_t i, const index_t j, const index_t k) const {
    return ((i + structure.gp) + ((j + structure.gp) + (k + structure.gp) * (structure.gy)) * (structure.gx));
}

void GridData::initData(const VectorFunction &initial_velocity, const Function &initial_pressure) {
    for (index_t z = 0; z < structure.nz; ++z) {
        for (index_t y = 0; y < structure.ny; ++y) {
            for (index_t x = 0; x < structure.nx; ++x) {
                Real px = real(x + structure.px) * structure.dx;
                Real py = real(y + structure.py) * structure.dy;
                Real pz = real(z + structure.pz) * structure.dz;

                U(x, y, z) = initial_velocity(px + structure.dx, py + structure.sdy, pz + structure.sdz)[0];
                V(x, y, z) = initial_velocity(px + structure.sdx, py + structure.dy, pz + structure.sdz)[1];
                W(x, y, z) = initial_velocity(px + structure.sdx, py + structure.sdy, pz + structure.dz)[2];
                P(x, y, z) = initial_pressure(px + structure.sdx, py + structure.sdy, pz + structure.sdz);
            }
        }
    }
}