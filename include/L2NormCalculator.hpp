#ifndef AEROHPC_A_L2NORM_CALCULATOR_H
#define AEROHPC_A_L2NORM_CALCULATOR_H

#include "StaggeredGrid.hpp"
#include <cmath>
#include <stdexcept>
#include <array>

// let's assume the spacing factor is h in all dimensions:
constexpr double h = 0.1;

// Define the exact solution functions for the u, v, and w components
class ExactSolution {
public:
    static double u(double x, double y, double z, double t) {
        return std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
    }

    static double v(double x, double y, double z, double t) {
        return std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
    }

    static double w(double x, double y, double z, double t) {
        return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
    }
};

class L2NormCalculator {
public:
    // Static method to compute the L2 norm of the difference between the grid values and exact solutions
    static double computeL2Norm(const /*StaggeredGrid<double, Addressing_T::A>*/& grid, double time) {
        double sum = 0.0;

        // Loop through the entire grid
        for (size_t k = 0; k < grid.getNz(); ++k) {
            for (size_t j = 0; j < grid.getNy(); ++j) {
                for (size_t i = 0; i < grid.getNx(); ++i) {

                    // Convert grid indices to real space coordinates
                    double x = static_cast<double>(i) * h;
                    double y = static_cast<double>(j) * h;
                    double z = static_cast<double>(k) * h;

                    // Calculate the exact solution for each component
                    double exactU = ExactSolution<double>::u(x, y, z, time);
                    double exactV = ExactSolution<double>::v(x, y, z, time);
                    double exactW = ExactSolution<double>::w(x, y, z, time);

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

        return std::sqrt(sum) / (h * h * h);
    }
};

#endif // AEROHPC_A_L2NORM_CALCULATOR_H
