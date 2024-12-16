#ifndef AEROHPC_A_TRAITS_H
#define AEROHPC_A_TRAITS_H

#include <functional>

/**
 * Typedef for real values
 */
typedef double Real;
#define real(val) static_cast<Real>(val)

/**
 * Typedef for array indexing
 */
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

/**
 * Typedef shortening Real Physical vector
 */
typedef std::array<Real, 3> Vector;

/**
 * Typedef shortening VolSpace dimensions
 */
typedef std::array<index_t, 3> Idx3;

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
 * Typedef shortening lambda definition of spatial function, time dependent
 */
typedef Real (*TFunction)(Real x, Real y, Real z, Real t);

/**
 * Typedef shortening lambda definition of spatial function, time independent
 */
typedef Real (*Function)(Real x, Real y, Real z);

/**
 * Typedef shortening lambda definition of vector spatial function
 */
typedef Vector (*VectorFunction)(Real x, Real y, Real z);


#endif //AEROHPC_A_TRAITS_H
