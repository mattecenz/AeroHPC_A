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


#include "mpi.h"

#if REAL_USE_FLOAT
#define Real_MPI MPI_FLOAT
#else
#define Real_MPI MPI_DOUBLE
#endif

#define NORTH_BUFFER_TAG 1
#define SOUTH_BUFFER_TAG 2
#define WEST_BUFFER_TAG 4
#define EAST_BUFFER_TAG 8
#define NORTH_EAST_BUFFER_TAG 16
#define NORTH_WEST_BUFFER_TAG 32
#define SOUTH_EAST_BUFFER_TAG 64
#define SOUTH_WEST_BUFFER_TAG 128

#define U_BUFFER_TAG 256
#define V_BUFFER_TAG 512
#define W_BUFFER_TAG 1024
#define P_BUFFER_TAG 2048


#endif //AEROHPC_A_TRAITS_H
