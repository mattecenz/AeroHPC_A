#include "L2NormCalculator.hpp"
#include <cmath> // For std::sqrt

Real ExactSolution::u(Real x, Real y, Real z, Real t) {
#ifdef TEST
    return std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
#else
    return 0;
#endif
}

Real ExactSolution::v(Real x, Real y, Real z, Real t) {
#ifdef TEST
    return std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
#else
    return 0;
#endif
}

Real ExactSolution::w(Real x, Real y, Real z, Real t) {
#ifdef TEST
    return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
#else
    return 0;
#endif
}


Real computeL2Norm(const GridData &grid, Real time) {
    Real sum = 0.0;
    Real sum_exact = 0.0;
    
    Real sdx = grid.structure.sdx;
    Real sdy = grid.structure.sdy;
    Real sdz = grid.structure.sdz;

    // Access the grid from the model
    // maybe change this later

    // Loop through the entire grid
    for (index_t i = 0; i < grid.structure.nx; ++i) {
        for (index_t j = 0; j < grid.structure.ny; ++j) {
            for (index_t k = 0; k < grid.structure.nz; ++k) {

                // Convert grid indices to real space coordinates
                Real x = real(i + grid.structure.px) * grid.structure.dx;
                Real y = real(j + grid.structure.py) * grid.structure.dy;
                Real z = real(k + grid.structure.pz) * grid.structure.dz;

                // Calculate the exact solution for each component
                Real exactU = ExactSolution::u(x + grid.structure.dx, y + sdy, z + sdz, time);
                Real exactV = ExactSolution::v(x + sdx, y + grid.structure.dy, z + sdz, time);
                Real exactW = ExactSolution::w(x + sdx, y + sdy, z + grid.structure.dz, time);

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

    return sum * (grid.structure.dx * grid.structure.dy * grid.structure.dz);
}
