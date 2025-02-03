#ifndef AEROHPC_A_TRAITS_H
#define AEROHPC_A_TRAITS_H

#include <functional>
#include <map>
#include <string>

#ifndef ForcingT
#define ForcingT 0
#endif

#ifndef DISABLE_PRESSURE
#define DISABLE_PRESSURE 0
#endif

#ifndef DEBUG_PRINT_BUFFERS
#define DEBUG_PRINT_BUFFERS 0
#endif


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
 * Typedef for solver result
 */
typedef std::map<std::string, Real> result_t;

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

/// Number of velocity components
constexpr unsigned int VELOCITY_COMPONENTS = 3;

/// Types of boundary condition
enum BOUNDARY_TYPE {
    DIRICHLET,
    NEUMANN
};

/// Problem unknowns
enum UNKNOWNS {
    VELOCITY = 1,
    PRESSURE = 2
};


#include "mpi.h"

inline int THIS_PROC_RANK = MPI_PROC_NULL;
inline int THIS_WORLD_SIZE = -1;

#define IS_MAIN_PROC (THIS_PROC_RANK == 0)

#if REAL_USE_FLOAT
#define Real_MPI MPI_FLOAT
#else
#define Real_MPI MPI_DOUBLE
#endif

enum BUFFER_TAGS {
    // Position buffer
    NORTH_BUFFER_TAG = 1,
    SOUTH_BUFFER_TAG = 2,
    WEST_BUFFER_TAG = 4,
    EAST_BUFFER_TAG = 8,
    NORTH_EAST_BUFFER_TAG = 16,
    NORTH_WEST_BUFFER_TAG = 32,
    SOUTH_EAST_BUFFER_TAG = 64,
    SOUTH_WEST_BUFFER_TAG = 128,
    // Component buffer
    U_BUFFER_TAG = 256,
    V_BUFFER_TAG = 512,
    W_BUFFER_TAG = 1024,
    P_BUFFER_TAG = 2048
};

#endif //AEROHPC_A_TRAITS_H
