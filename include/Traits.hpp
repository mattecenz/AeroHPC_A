#ifndef AEROHPC_A_TRAITS_H
#define AEROHPC_A_TRAITS_H

#include <functional>

#if REAL_USE_FLOAT
typedef float Real;
#define Real_Dataname "float"
#else
typedef double Real;
#define Real_Dataname "double"
#endif

#define real(val) static_cast<Real>(val)
#define real_p(pntr) static_cast<Real*>(pntr)

/**
 * Typedef for array indexing
 */
typedef long index_t;

/**
 * Typedef shortening Real Physical vector
 */
typedef std::array<Real, 3> Vector;

/**
 * Typedef shortening VolSpace dimensions
 */
typedef std::array<index_t, 3> Idx3;

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
