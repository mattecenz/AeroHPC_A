#include "L2NormCalculator.hpp"
#include <cmath> // For std::sqrt

template<Addressing_T A>
Real computeL2Norm(const Model<A> &model, Real time) {
    Real sum = 0.0;

    Real sdx = model.sdx;
    Real sdy = model.sdy;
    Real sdz = model.sdz;

    // Access the grid from the model
    // maybe change this later
    const auto &grid = model.grid; // Access grid from model

    // Loop through the entire grid
    for (index_t i = 0; i < grid.nx; ++i) {
        for (index_t j = 0; j < grid.ny; ++j) {
            for (index_t k = 0; k < grid.nz; ++k) {

                // Convert grid indices to real space coordinates
                Real x = real(i) * model.dx;
                Real y = real(j) * model.dy;
                Real z = real(k) * model.dz;

                // Calculate the exact solution for each component
                Real exactU = ExactSolution<A>::u(x + model.dx, y + sdy, z + sdz, time);
                Real exactV = ExactSolution<A>::v(x + sdx, y + model.dy, z + sdz, time);
                Real exactW = ExactSolution<A>::w(x + sdx, y + sdy, z + model.dz, time);

                // Access the computed grid components
                Real gridU = grid(Component::U, i, j, k);
                Real gridV = grid(Component::V, i, j, k);
                Real gridW = grid(Component::W, i, j, k);

                // Calculate the differences
                Real diffU = gridU - exactU;
                Real diffV = gridV - exactV;
                Real diffW = gridW - exactW;

                // Add the squares of the differences to sum
                sum += (diffU * diffU) + (diffV * diffV) + (diffW * diffW);
            }
        }
    }

    return std::sqrt(sum * (model.dx * model.dy * model.dz));
}

// Explicit instantiation for the Addressing_T 
template Real computeL2Norm<STANDARD>(const Model<STANDARD> &model, Real time);
