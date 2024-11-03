#include "Grid.hpp"

template<>
Real &Grid<STANDARD>::operator()(const Component c, const index_t i, const index_t j, const index_t k) {
    return _entries[c][((i + 1) + ((j + 1) + (k + 1) * (ny + 2 * 1)) * (ny + 2 * 1))];
}

template<>
Real Grid<STANDARD>::operator()(const Component c, const index_t i, const index_t j, const index_t k) const  {
    return _entries[c][((i + 1) + ((j + 1) + (k + 1) * (ny + 2 * 1)) * (ny + 2 * 1))];
}

template<>
void Grid<STANDARD>::initGrid(const VectorFunction &initial_velocity, const Function &initial_pressure) {
    for (index_t z = 0; z < nz; ++z) {
        for (index_t y = 0; y < ny; ++y) {
            for (index_t x = 0; x < nx; ++x) {
                Real px = real(x) * dx;
                Real py = real(y) * dy;
                Real pz = real(z) * dz;

                operator()(U, x, y, z) = initial_velocity(px + dx, py + sdy, pz + sdz)[0];
                operator()(V, x, y, z) = initial_velocity(px + sdx, py + dy, pz + sdz)[1];
                operator()(W, x, y, z) = initial_velocity(px + sdx, py + sdy, pz + dz)[2];
                operator()(P, x, y, z) = initial_pressure(px + sdx, py + sdy, pz + sdz);
            }
        }
    }
}