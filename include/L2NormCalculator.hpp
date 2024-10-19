#ifndef AEROHPC_A_L2NORM_CALCULATOR_H
#define AEROHPC_A_L2NORM_CALCULATOR_H

#include <cmath>
#include <stdexcept>
#include <array>
#include "StaggeredGrid.hpp"
#include "utils.hpp" 

// Define the exact solution functions for the u, v, and w components
template<typename T>
class ExactSolution {
public:
    static T u(T x, T y, T z, T t) {
        return std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
    }

    static T v(T x, T y, T z, T t) {
        return std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
    }

    static T w(T x, T y, T z, T t) {
        return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
    }
};


template<typename T>
class L2NormCalculator {
public:
    static typename Traits<T>::Real computeL2Norm(const utils::StaggeredGridType<T>& grid, T time) {
        typename Traits<T>::Real sum = 0.0;
        
        T h = 0.1;

        for (size_t k = 0; k < grid.getNz(); ++k) {
            for (size_t j = 0; j < grid.getNy(); ++j) {
                for (size_t i = 0; i < grid.getNx(); ++i) {
                    T x = static_cast<T>(i) * h;
                    T y = static_cast<T>(j) * h;
                    T z = static_cast<T>(k) * h;

                    // Calculate exact solution
                    T exactU = ExactSolution<T>::u(x, y, z, time);
                    T exactV = ExactSolution<T>::v(x, y, z, time);
                    T exactW = ExactSolution<T>::w(x, y, z, time);

                    // Access grid components
                    T gridU = grid(Component::U, i, j, k);
                    T gridV = grid(Component::V, i, j, k);
                    T gridW = grid(Component::W, i, j, k);

                    // Compute differences
                    T diffU = gridU - exactU;
                    T diffV = gridV - exactV;
                    T diffW = gridW - exactW;

                    // Accumulate the sum of squares
                    sum += (diffU * diffU) + (diffV * diffV) + (diffW * diffW);
                }
            }
        }

        // Return the L2 norm
        return std::sqrt(sum) / (h * h * h);
    }
};

#endif // AEROHPC_A_L2NORM_CALCULATOR_H
