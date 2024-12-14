#ifndef AEROHPC_A_MPITRAITS_HPP
#define AEROHPC_A_MPITRAITS_HPP

#include "mpi.h"

#if Real == float
#define Real_MPI MPI_FLOAT
#elif Real == double
#define Real_MPI MPI_DOUBLE
#endif

#define NORTH_BUFFER_TAG 0
#define SOUTH_BUFFER_TAG 1
#define WEST_BUFFER_TAG 2
#define EAST_BUFFER_TAG 3
#define NORTH_EAST_BUFFER_TAG 4
#define NORTH_WEST_BUFFER_TAG 5
#define SOUTH_EAST_BUFFER_TAG 6
#define SOUTH_WEST_BUFFER_TAG 7

#endif //AEROHPC_A_MPITRAITS_HPP
