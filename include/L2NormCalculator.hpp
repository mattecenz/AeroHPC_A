#ifndef AEROHPC_A_L2NORM_CALCULATOR_H
#define AEROHPC_A_L2NORM_CALCULATOR_H

#include "Traits.hpp"

// Define the exact solution functions for the u, v, and w components
class ExactSolution {
public:
    static Real u(Real x, Real y, Real z, Real t);

    static Real v(Real x, Real y, Real z, Real t);

    static Real w(Real x, Real y, Real z, Real t);
};

Real computeL2Norm(const Real *grid, Real time);

#endif // AEROHPC_A_L2NORM_CALCULATOR_H
