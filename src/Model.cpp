
#include "Model.hpp"

/**
 * Initialize the space grid with standard addressing
 */
template<>
void Model<Addressing_T::STANDARD>::initGrid(const Function& initial_velocity, const Function& initial_pressure) {
    Real sdx = _spacing[0]/2;
    Real sdy = _spacing[1]/2;
    Real sdz = _spacing[2]/2;

    for (size_t z = 0; z < _grid._nz; ++z) {
        for (size_t y = 0; y < _grid._ny; ++y) {
            for (size_t x = 0; x < _grid._nx; ++x) {
                _grid(U, x, y, z) = initial_velocity(x + sdx, y, z);
                _grid(V, x, y, z) = initial_velocity(x, y + sdy, z);
                _grid(W, x, y, z) = initial_velocity(x, y, z + sdz);
                _grid(P, x, y, z) = initial_pressure(x, y, z);
            }
        }
    }
}