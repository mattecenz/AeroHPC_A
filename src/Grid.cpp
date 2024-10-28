#include "Grid.hpp"

template<>
Real &Grid<STANDARD>::operator()(Component c, index_t i, index_t j, index_t k) {
    return _entries[c + ((i + _gp) + ((j + _gp) + (k + _gp) * _grid_nodes[1]) * _grid_nodes[0]) * N_COMPONENTS];
}

template<>
const Real &Grid<STANDARD>::operator()(Component c, index_t i, index_t j, index_t k) const {
    return _entries[c + ((i + _gp) + ((j + _gp) + (k + _gp) * _grid_nodes[1]) * _grid_nodes[0]) * N_COMPONENTS];
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