#ifndef AEROHPC_A_MPITRAITS_HPP
#define AEROHPC_A_MPITRAITS_HPP

#include "mpi.h"

#if Real == float
#define MPI_Real MPI_FLOAT
#elif Real == double
#define MPI_Real MPI_DOUBLE
#endif

#define NORTH_BUFFER_TAG 0
#define SOUTH_BUFFER_TAG 1
#define WEST_BUFFER_TAG 2
#define EAST_BUFFER_TAG 3

#endif //AEROHPC_A_MPITRAITS_HPP
