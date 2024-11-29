#ifndef AEROHPC_A_L2NORM_CALCULATOR_H
#define AEROHPC_A_L2NORM_CALCULATOR_H

#include "GridData.hpp"
#include <cmath>
#include <stdexcept>
#include <array>

// let's assume the spacing factor is h in all dimensions:
// constexpr Real h = 0.1;

// Define the exact solution functions for the u, v, and w components
class ExactSolution {
public:
    static Real u(Real x, Real y, Real z, Real t) {
        return std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
    }

    static Real v(Real x, Real y, Real z, Real t) {
        return std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
    }

    static Real w(Real x, Real y, Real z, Real t) {
        return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
    }
};

Real computeL2Norm(const GridData &grid, Real time);

#endif // AEROHPC_A_L2NORM_CALCULATOR_H
