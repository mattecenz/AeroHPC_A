#ifndef AEROHPC_A_L2NORM_CALCULATOR_H
#define AEROHPC_A_L2NORM_CALCULATOR_H

#include "StaggeredGrid.hpp"
#include "Model.hpp"
#include <cmath>
#include <stdexcept>
#include <array>

// let's assume the spacing factor is h in all dimensions:
// constexpr Real h = 0.1;

// Define the exact solution functions for the u, v, and w components
template <Addressing_T A>
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

template <Addressing_T A>
Real computeL2Norm(const Model<A> &model, Real time);

#endif // AEROHPC_A_L2NORM_CALCULATOR_H
