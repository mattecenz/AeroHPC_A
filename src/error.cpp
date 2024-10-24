#include "L2NormCalculator.hpp"
#include <cmath> // For std::sqrt

template <Addressing_T A>
double computeL2Norm(const Model<A> &model, double time) {
    double sum = 0.0;

    // Access the grid from the model
    // maybe change this later
    const auto &grid = model.grid; // Access grid from model

    // Loop through the entire grid
    for (size_t k = 0; k < grid.nx; ++k) {
        for (size_t j = 0; j < grid.ny; ++j) {
            for (size_t i = 0; i < grid.nz; ++i) {

                // Convert grid indices to real space coordinates
                double x = static_cast<double>(i) * model.dx;
                double y = static_cast<double>(j) * model.dy;
                double z = static_cast<double>(k) * model.dz;

                // Calculate the exact solution for each component
                double exactU = ExactSolution<A>::u(x, y, z, time);
                double exactV = ExactSolution<A>::v(x, y, z, time);
                double exactW = ExactSolution<A>::w(x, y, z, time);

                // Access the computed grid components
                double gridU = grid(Component::U, i, j, k);
                double gridV = grid(Component::V, i, j, k);
                double gridW = grid(Component::W, i, j, k);

                // Calculate the differences
                double diffU = gridU - exactU;
                double diffV = gridV - exactV;
                double diffW = gridW - exactW;

                // Add the squares of the differences to sum
                sum += (diffU * diffU) + (diffV * diffV) + (diffW * diffW);
            }
        }
    }

    return std::sqrt(sum) / (model.dx * model.dy * model.dz);
}

// Explicit instantiation for the Addressing_T 
template double computeL2Norm<STANDARD>(const Model<STANDARD> &model, double time);
