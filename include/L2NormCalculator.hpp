#ifndef AEROHPC_A_L2NORM_CALCULATOR_H
#define AEROHPC_A_L2NORM_CALCULATOR_H

#include <StaggeredGrid.hpp>
#include <cmath>
#include <stdexcept>
#include <array>

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


template<typename T, Addressing_T A>
class L2NormCalculator {
public:
    // Constructor
    L2NormCalculator(const StaggeredGrid<T, Addressing_T::A>& grid, T time)
        : _grid(grid), _time(time) {}

    // Function to compute the L2 norm of the difference between the grid values and exact solutions
    T computeL2Norm() const {
        T sum = 0.0;

        // Loop through the entire grid
        for (size_t k = 0; k < _grid.getNz(); ++k) {
            for (size_t j = 0; j < _grid.getNy(); ++j) {
                for (size_t i = 0; i < _grid.getNx(); ++i) {

                    // let's assume the spacing factor is h in all dimensions:
                    // ATENTION //////////////////
                    // change this:
                    T h = 0.1; 

                    // Convert grid indices to real space coordinates
                    T x = static_cast<T>(i) * h;
                    T y = static_cast<T>(j) * h;
                    T z = static_cast<T>(k) * h;

                    // Calculate the exact solution for each component
                    T exactU = ExactSolution<T>::u(x, y, z, _time);
                    T exactV = ExactSolution<T>::v(x, y, z, _time);
                    T exactW = ExactSolution<T>::w(x, y, z, _time);

                    // Access the computed grid components
                    T gridU = _grid(Component::U, i, j, k);
                    T gridV = _grid(Component::V, i, j, k);
                    T gridW = _grid(Component::W, i, j, k);

                    // Calculate the differences
                    T diffU = gridU - exactU;
                    T diffV = gridV - exactV;
                    T diffW = gridW - exactW;

                    // Add the squares of the differences to sum
                    sum += (diffU * diffU) + (diffV * diffV) + (diffW * diffW);


                }
            }
        }

        return std::sqrt(sum) / (h * h * h);
    }

private:
    const StaggeredGrid<T, Addressing_T::STANDARD>& _grid; // Reference to the grid
    T _time;                          // Current time value
};

#endif // AEROHPC_A_L2NORM_CALCULATOR_H
