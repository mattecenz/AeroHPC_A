#ifndef AEROHPC_A_TRAITS_H
#define AEROHPC_A_TRAITS_H

#include <functional>

/**
 * Typedef for real values
 */
typedef float Real;

typedef long index_t;

/**
 * Enum for easy modification of the class addressing
 */
enum Addressing_T {
    /**
     * Multidimensional grid flattened into a single dimension array composed by cells of 4 values:
     * [u000, v000, w000, p000, u100, v100, w100, p100,  ... ]
     */
    STANDARD = 0
};

typedef std::array<Real, 3> Vector;

/**
 * Define mathematical methods for vectors
 */
#pragma inline

Vector operator*(Real b, const Vector &a);

#pragma inline

Vector operator*(const Vector &a);

#pragma inline

Real operator*(const Vector &a, const Vector &b);

#pragma inline

Vector operator+(const Vector &a, const Vector &b);

#pragma inline

Vector operator-(const Vector &a);


/**
 * Typedef shortening lambda definition of spatial function
 */
typedef std::function<Real(Real x, Real y, Real z)> Function;

/**
 * Typedef shortening lambda definition of vector spatial function
 */
typedef std::function<Vector(Real x, Real y, Real z)> VectorFunction;

#endif //AEROHPC_A_TRAITS_H
