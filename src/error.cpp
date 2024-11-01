#include "L2NormCalculator.hpp"
#include <cmath> // For std::sqrt

Real computeL2Norm(const Grid &grid, Real time) {
    Real sum = 0.0;
    Real sum_exact = 0.0;
    
    Real sdx = grid.sdx;
    Real sdy = grid.sdy;
    Real sdz = grid.sdz;

    // Access the grid from the model
    // maybe change this later

    // Loop through the entire grid
    for (index_t i = 0; i < grid.nx; ++i) {
        for (index_t j = 0; j < grid.ny; ++j) {
            for (index_t k = 0; k < grid.nz; ++k) {

                // Convert grid indices to real space coordinates
                Real x = real(i) * grid.dx;
                Real y = real(j) * grid.dy;
                Real z = real(k) * grid.dz;

                // Calculate the exact solution for each component
                Real exactU = ExactSolution::u(x + grid.dx, y + sdy, z + sdz, time);
                Real exactV = ExactSolution::v(x + sdx, y + grid.dy, z + sdz, time);
                Real exactW = ExactSolution::w(x + sdx, y + sdy, z + grid.dz, time);

                // Access the computed grid components
                Real gridU = grid.U(i, j, k);
                Real gridV = grid.V(i, j, k);
                Real gridW = grid.W(i, j, k);

                // Calculate the differences
                Real diffU = gridU - exactU;
                Real diffV = gridV - exactV;
                Real diffW = gridW - exactW;

                // Add the squares of the differences to sum
                sum += (diffU * diffU) + (diffV * diffV) + (diffW * diffW);
                sum_exact += (exactU * exactU) + (exactV * exactV) + (exactW * exactW);
            }
        }
    }

    return std::sqrt(sum * (grid.dx * grid.dy * grid.dz));///std::sqrt(sum_exact * (grid.dx * grid.dy * grid.dz));
}
