/**
 * This file contains functions that compile the BC collection
 */

#include "C2Decomp.hpp"

#include "Traits.hpp"
#include "Condition.hpp"
#include "Boundaries.hpp"

#include "MPITraits.hpp"
#include "MPICondition.hpp"
#include "MPIBoundaries.hpp"

// Couple of macro for shortening the code
#define getStaggeredSpacing(grid, x, y, z) \
    const Real x = grid.structure.sdx; \
    const Real y = grid.structure.sdy; \
    const Real z = grid.structure.sdz

#define getExactFunctions(functions, u, v, w) \
    const TFunction u = functions[0]; \
    const TFunction v = functions[1]; \
    const TFunction w = functions[2]

#define cast_int(v) static_cast<int>(v)

#define NORTH_FACE_ID 0
#define SOUTH_FACE_ID 1
#define EAST_FACE_ID 2
#define WEST_FACE_ID 3
#define FRONT_FACE_ID 4
#define BACK_FACE_ID 5

typedef std::vector<TFunction> boundaryFaceFunctions;
typedef std::array<boundaryFaceFunctions, 6> boundaryDomainFunctions;

/// NORTH //////////////////////////////////////////////////////////////////////////////////////////////
namespace north {
    // This lambda defines how the #functions have to be applied to the #grid
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const boundaryFaceFunctions &functions) {
        // The upper boundary is at the max Y
        const index_t j = grid.structure.ny;

        // PERIODIC CONDITION
        if (functions.empty()) {
            const index_t periodic_j = 0;

            // apply on face
            for (index_t k = 0; k < grid.structure.nz; k++)
                for (index_t i = 0; i < grid.structure.nx; i++) {
                    grid.U(i, j, k) = grid.U(i, periodic_j, k);
                    grid.V(i, j, k) = grid.V(i, periodic_j, k);
                    grid.W(i, j, k) = grid.W(i, periodic_j, k);

                    grid.P(i, j, k) = grid.P(i, j - 1, k);
                }
        }
        // DIRICHLET CONDITION
        else {
            // Use macro to get some variables
            getStaggeredSpacing(grid, sdx, sdy, sdz);
            getExactFunctions(functions, eU, eV, eW);

            const Real y = real(j + grid.structure.py) * grid.structure.dy;

            // apply on face
            for (index_t k = 0; k < grid.structure.nz; k++)
                for (index_t i = 0; i < grid.structure.nx; i++) {
                    Real x = real(i + grid.structure.px) * grid.structure.dx;
                    Real z = real(k + grid.structure.pz) * grid.structure.dz;

                    // On y = phy_dim for domain point we have exact for V
                    grid.V(i, j - 1, k) = eV(x + sdx, y, z + sdz, currentTime);
                    // For ghost points we have useless V, other approximate

                    grid.U(i, j, k) = 2 * eU(x + grid.structure.dx, y, z + sdz, currentTime)
                                      - grid.U(i, j - 1, k);
                    grid.V(i, j, k) = 0;
                    grid.W(i, j, k) = 2 * eW(x + sdx, y, z + grid.structure.dz, currentTime)
                                      - grid.W(i, j - 1, k);

                    // P derivative has to be 0 in y direction, so the value should be the same
                    grid.P(i, j, k) = grid.P(i, j - 1, k);
                }
        }
    };
    // This lambda defines how the outgoing communication buffer has to be initialized
    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t j = grid.structure.ny - 1;
        for (index_t k = -1; k <= grid.structure.nz; k++) {
            memcpy(&bufferOut.U(0, 0, k), &grid.U(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.V(0, 0, k), &grid.V(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.W(0, 0, k), &grid.W(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.P(0, 0, k), &grid.P(-1, j, k), sizeof(Real) * grid.structure.gx);
        }
    };

    // This lambda defines which data has to be shared with the given neighbour (#neigh_rank)
    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        // This proc will send his outgoing buffer with the tag #NORTH_BUFFER_TAG
        // (means that the buffer is the top layer of this domain)
        MPI_Isend(bufferOut.data, cast_int(bufferOut.node_dim) * 4, Real_MPI,
                  neigh_rank, NORTH_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        // This proc will receive into ingoing buffer data with tag #SOUTH_BUFFER_TAG
        // (means that the buffer is the top layer of neighbour domain)
        MPI_Irecv(bufferIn.data, cast_int(bufferIn.node_dim) * 4, Real_MPI,
                  neigh_rank, SOUTH_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    // This lambda defines how the ingoing communication buffer has to be copied into the local model grid
    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        // Since communication are asynchronous I wait for the ingoing buffer to be completely received
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);

        // Copy the buffer into the upper ghost point layer
        const index_t j = grid.structure.ny;
        for (index_t k = -1; k <= grid.structure.nz; k++) {
            memcpy(&grid.U(-1, j, k), &buffer.U(0, 0, k), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.V(-1, j, k), &buffer.V(0, 0, k), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.W(-1, j, k), &buffer.W(0, 0, k), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.P(-1, j, k), &buffer.P(0, 0, k), sizeof(Real) * grid.structure.gx);
        }
    };
}

/// SOUTH //////////////////////////////////////////////////////////////////////////////////////////////
namespace south {
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const boundaryFaceFunctions &functions) {
        const index_t j = -1;

        if (functions.empty()) {
            const index_t periodic_j = grid.structure.ny - 1;

            for (index_t k = 0; k < grid.structure.nz; k++)
                for (index_t i = 0; i < grid.structure.nx; i++) {
                    grid.U(i, j, k) = grid.U(i, periodic_j, k);
                    grid.V(i, j, k) = grid.V(i, periodic_j, k);
                    grid.W(i, j, k) = grid.W(i, periodic_j, k);
                    grid.P(i, j - 1, k) = grid.P(i, j, k);
                }
        } else {
            getStaggeredSpacing(grid, sdx, sdy, sdz);
            getExactFunctions(functions, eU, eV, eW);

            // South face has y=0
            const Real y = real((j+1) + grid.structure.py) * grid.structure.dy;

            for (index_t k = 0; k < grid.structure.nz; k++)
                for (index_t i = 0; i < grid.structure.nx; i++) {
                    Real x = real(i + grid.structure.px) * grid.structure.dx;
                    Real z = real(k + grid.structure.pz) * grid.structure.dz;

                    // On y = 0 for ghost point we hae exact for V, other approximate
                    grid.U(i, j, k) = 2 * eU(x + grid.structure.dx, y, z + sdz, currentTime)
                                      - grid.U(i, j + 1, k);
                    grid.V(i, j, k) = eV(x + sdx, y, z + sdz, currentTime);
                    grid.W(i, j, k) = 2 * eW(x + sdx, y, z + grid.structure.dz, currentTime)
                                      - grid.W(i, j + 1, k);

                    grid.P(i, j, k) = grid.P(i, j + 1, k);
                }
        }
    };

    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t j = 0;
        for (index_t k = -1; k <= grid.structure.nz; k++) {
            memcpy(&bufferOut.U(0, 0, k), &grid.U(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.V(0, 0, k), &grid.V(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.W(0, 0, k), &grid.W(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.P(0, 0, k), &grid.P(-1, j, k), sizeof(Real) * grid.structure.gx);
        }
    };

    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        MPI_Isend(bufferOut.data, cast_int(bufferOut.node_dim) * 4, Real_MPI,
                  neigh_rank, SOUTH_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        MPI_Irecv(bufferIn.data, cast_int(bufferIn.node_dim) * 4, Real_MPI,
                  neigh_rank, NORTH_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);

        // Copy the buffer into the ghost point layer
        const index_t j = -1;
        for (index_t k = -1; k <= grid.structure.nz; k++) {
            memcpy(&grid.U(-1, j, k), &buffer.U(0, 0, k), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.V(-1, j, k), &buffer.V(0, 0, k), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.W(-1, j, k), &buffer.W(0, 0, k), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.P(-1, j, k), &buffer.P(0, 0, k), sizeof(Real) * grid.structure.gx);
        }
    };
}

/// EAST ///////////////////////////////////////////////////////////////////////////////////////////////
namespace east {
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const boundaryFaceFunctions &functions) {
        const index_t k = grid.structure.nz;

        if (functions.empty()) {
            const index_t periodic_k = 0;

            for (index_t j = 0; j < grid.structure.ny; j++)
                for (index_t i = 0; i < grid.structure.nx; i++) {
                    grid.U(i, j, k) = grid.U(i, j, periodic_k);
                    grid.V(i, j, k) = grid.V(i, j, periodic_k);
                    grid.W(i, j, k) = grid.W(i, j, periodic_k);
                    grid.P(i, j, k) = grid.P(i, j, k - 1);
                }
        } else {
            getStaggeredSpacing(grid, sdx, sdy, sdz);
            getExactFunctions(functions, eU, eV, eW);

            const Real z = real(k + grid.structure.pz) * grid.structure.dz;

            for (index_t j = 0; j < grid.structure.ny; j++)
                for (index_t i = 0; i < grid.structure.nx; i++) {
                    Real x = real(i + grid.structure.px) * grid.structure.dx;
                    Real y = real(j + grid.structure.py) * grid.structure.dy;

                    // On z = phy_dim for domain point we have exact for W
                    grid.W(i, j, k - 1) = eW(x + sdx, y + sdy, z, currentTime);
                    // For ghost points we gave useless W, other interpolate
                    grid.U(i, j, k) = 2 * eU(x + grid.structure.dx, y + sdy, z, currentTime)
                                      - grid.U(i, j, k - 1);
                    grid.V(i, j, k) = 2 * eV(x + sdx, y + grid.structure.dy, z, currentTime)
                                      - grid.V(i, j, k - 1);
                    grid.W(i, j, k) = 0;

                    grid.P(i, j, k) = grid.P(i, j, k - 1);
                }
        }
    };

    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t k = grid.structure.nz - 1;

        for (index_t j = -1; j <= grid.structure.ny; j++) {
            memcpy(&bufferOut.U(0, j, 0), &grid.U(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.V(0, j, 0), &grid.V(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.W(0, j, 0), &grid.W(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.P(0, j, 0), &grid.P(-1, j, k), sizeof(Real) * grid.structure.gx);
        }
    };

    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        MPI_Isend(bufferOut.data, cast_int(bufferOut.node_dim) * 4, Real_MPI, neigh_rank,
                  EAST_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        MPI_Irecv(bufferIn.data, cast_int(bufferIn.node_dim) * 4, Real_MPI, neigh_rank,
                  WEST_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);

        // Copy the buffer into the ghost point layer
        const index_t k = grid.structure.nz;
        for (index_t j = -1; j <= grid.structure.ny; j++) {
            memcpy(&grid.U(-1, j, k), &buffer.U(0, j, 0), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.V(-1, j, k), &buffer.V(0, j, 0), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.W(-1, j, k), &buffer.W(0, j, 0), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.P(-1, j, k), &buffer.P(0, j, 0), sizeof(Real) * grid.structure.gx);
        }
    };
}

/// WEST ///////////////////////////////////////////////////////////////////////////////////////////////
namespace west {
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const boundaryFaceFunctions &functions) {
        const index_t k = -1;

        if (functions.empty()) {
            const index_t periodic_k = grid.structure.nz - 1;

            for (index_t j = 0; j < grid.structure.ny; j++)
                for (index_t i = 0; i < grid.structure.nx; i++) {
                    grid.U(i, j, k) = grid.U(i, j, periodic_k);
                    grid.V(i, j, k) = grid.V(i, j, periodic_k);
                    grid.W(i, j, k) = grid.W(i, j, periodic_k);
                    grid.P(i, j, k) = grid.P(i, j, k + 1);
                }

        } else {
            getStaggeredSpacing(grid, sdx, sdy, sdz);
            getExactFunctions(functions, eU, eV, eW);

            const Real z = real((k+1) + grid.structure.pz) * grid.structure.dz;

            for (index_t j = 0; j < grid.structure.ny; j++)
                for (index_t i = 0; i < grid.structure.nx; i++) {
                    Real x = real(i + grid.structure.px) * grid.structure.dx;
                    Real y = real(j + grid.structure.py) * grid.structure.dy;

                    // On z = 0 for ghost point we have exact for W, other approximate
                    grid.U(i, j, k) = 2 * eU(x + grid.structure.dx, y + sdy, z, currentTime)
                                      - grid.U(i, j, k + 1);
                    grid.V(i, j, k) = 2 * eV(x + sdx, y + grid.structure.dy, z, currentTime)
                                      - grid.V(i, j, k + 1);
                    grid.W(i, j, k) = eW(x + sdx, y + sdy, z, currentTime);

                    grid.P(i, j, k) = grid.P(i, j, k + 1);
                }
        }
    };

    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t k = 0;
        for (index_t j = -1; j <= grid.structure.ny; j++) {
            memcpy(&bufferOut.U(0, j, 0), &grid.U(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.V(0, j, 0), &grid.V(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.W(0, j, 0), &grid.W(-1, j, k), sizeof(Real) * grid.structure.gx);
            memcpy(&bufferOut.P(0, j, 0), &grid.P(-1, j, k), sizeof(Real) * grid.structure.gx);
        }
    };

    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        MPI_Isend(bufferOut.data, cast_int(bufferOut.node_dim) * 4, Real_MPI, neigh_rank,
                  WEST_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        MPI_Irecv(bufferIn.data, cast_int(bufferIn.node_dim) * 4, Real_MPI, neigh_rank,
                  EAST_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);

        // Copy the buffer into the ghost point layer
        const index_t k = -1;
        for (index_t j = -1; j <= grid.structure.ny; j++) {
            memcpy(&grid.U(-1, j, k), &buffer.U(0, j, 0), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.V(-1, j, k), &buffer.V(0, j, 0), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.W(-1, j, k), &buffer.W(0, j, 0), sizeof(Real) * grid.structure.gx);
            memcpy(&grid.P(-1, j, k), &buffer.P(0, j, 0), sizeof(Real) * grid.structure.gx);
        }
    };
}

/// FRONT //////////////////////////////////////////////////////////////////////////////////////////////
namespace front {
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const boundaryFaceFunctions &functions) {
        const index_t i = -1;

        if (functions.empty()) {
            const index_t periodic_i = grid.structure.nx - 1;

            for (index_t k = 0; k < grid.structure.nz; k++)
                for (index_t j = 0; j < grid.structure.ny; j++) {
                    grid.U(i,j,k) = grid.U(periodic_i, j, k);
                    grid.V(i,j,k) = grid.V(periodic_i, j, k);
                    grid.W(i,j,k) = grid.W(periodic_i, j, k);
                    grid.P(i, j, k) = grid.P(i, j+1, k);
                }
        }
        else{
            getStaggeredSpacing(grid, sdx, sdy, sdz);
            getExactFunctions(functions, eU, eV, eW);

            const Real x = real((i+1) + grid.structure.px) * grid.structure.dx;

            for (index_t k = 0; k < grid.structure.nz; k++)
                for (index_t j = 0; j < grid.structure.ny; j++) {
                    Real y = real(j + grid.structure.py) * grid.structure.dy;
                    Real z = real(k + grid.structure.pz) * grid.structure.dz;

                    // On x = 0 for ghost point we have exact for U, other approximate
                    grid.U(i, j, k) = eU(x, y + sdy, z + sdz, currentTime);
                    grid.V(i, j, k) = 2 * eV(x, y + grid.structure.dy, z + sdz, currentTime)
                                          - grid.V(i, j+1, k);
                    grid.W(i, j, k) = 2 * eW(x, y + sdy, z + grid.structure.dz, currentTime)
                                          - grid.W(i, j+1, k);

                    grid.P(i, j, k) = grid.P(i, j+1, k);
                }
        }
    };
}

/// BACK ///////////////////////////////////////////////////////////////////////////////////////////////
namespace back {
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const boundaryFaceFunctions &functions) {

        const index_t i = grid.structure.nx;

        if (functions.empty()) {
            const index_t periodic_i = 0;

            for (index_t k = 0; k < grid.structure.nz; k++)
                for (index_t j = 0; j < grid.structure.ny; j++) {
                    grid.U(i, j, k) = grid.U(periodic_i, j, k);
                    grid.V(i, j, k) = grid.V(periodic_i, j, k);
                    grid.W(i, j, k) = grid.W(periodic_i, j, k);

                    grid.P(i - 1, j, k) = grid.P(i, j, k);
                }
        }
        else {
            getStaggeredSpacing(grid, sdx, sdy, sdz);
            getExactFunctions(functions, eU, eV, eW);

            const Real x = real(i + grid.structure.px) * grid.structure.dx;

            for (index_t k = 0; k < grid.structure.nz; k++)
                for (index_t j = 0; j < grid.structure.ny; j++) {
                    Real y = real(j + grid.structure.py) * grid.structure.dy;
                    Real z = real(k + grid.structure.pz) * grid.structure.dz;

                    // On x = phy_dim for domain point we have exact for U
                    grid.U(i - 1, j, k) = eU(x, y + sdy, z + sdz, currentTime);
                    // For ghost point we have useless U, other approximate
                    grid.U(i, j, k) = 0;
                    grid.V(i, j, k) = 2 * eV(x, y + grid.structure.dy, z + sdz, currentTime)
                                      - grid.V(i - 1, j, k);
                    grid.W(i, j, k) = 2 * eW(x, y + sdy, z + grid.structure.dz, currentTime)
                                      - grid.W(i - 1, j, k);

                    grid.P(i - 1, j, k) = grid.P(i, j, k);
                }
        }
    };
}

/// NORTH EAST /////////////////////////////////////////////////////////////////////////////////////////
namespace nhet {
    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t j = grid.structure.ny - 1;
        const index_t k = grid.structure.nz - 1;
        memcpy(&bufferOut.U(0, 0, 0), &grid.U(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.V(0, 0, 0), &grid.V(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.W(0, 0, 0), &grid.W(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.P(0, 0, 0), &grid.P(-1, j, k), sizeof(Real) * grid.structure.gx);
    };

    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        MPI_Isend(bufferOut.data, cast_int(bufferOut.node_dim) * 4, Real_MPI,
                  neigh_rank, NORTH_EAST_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        MPI_Irecv(bufferIn.data, cast_int(bufferIn.node_dim) * 4, Real_MPI,
                  neigh_rank, SOUTH_WEST_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);
        const index_t j = grid.structure.ny;
        const index_t k = grid.structure.nz;
        memcpy(&grid.U(-1, j, k), &buffer.U(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.V(-1, j, k), &buffer.V(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.W(-1, j, k), &buffer.W(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.P(-1, j, k), &buffer.P(0, 0, 0), sizeof(Real) * grid.structure.gx);
    };
}

/// NORTH WEST /////////////////////////////////////////////////////////////////////////////////////////
namespace nhwt {
    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t j = grid.structure.ny - 1;
        const index_t k = 0;
        memcpy(&bufferOut.U(0, 0, 0), &grid.U(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.V(0, 0, 0), &grid.V(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.W(0, 0, 0), &grid.W(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.P(0, 0, 0), &grid.P(-1, j, k), sizeof(Real) * grid.structure.gx);
    };

    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        MPI_Isend(bufferOut.data, cast_int(bufferOut.node_dim) * 4, Real_MPI,
                  neigh_rank, NORTH_WEST_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        MPI_Irecv(bufferIn.data, cast_int(bufferIn.node_dim) * 4, Real_MPI,
                  neigh_rank, SOUTH_EAST_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);
        const index_t j = grid.structure.ny;
        const index_t k = -1;
        memcpy(&grid.U(-1, j, k), &buffer.U(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.V(-1, j, k), &buffer.V(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.W(-1, j, k), &buffer.W(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.P(-1, j, k), &buffer.P(0, 0, 0), sizeof(Real) * grid.structure.gx);
    };
}

/// SOUTH EAST /////////////////////////////////////////////////////////////////////////////////////////
namespace shet {
    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t j = 0;
        const index_t k = grid.structure.nz - 1;
        memcpy(&bufferOut.U(0, 0, 0), &grid.U(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.V(0, 0, 0), &grid.V(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.W(0, 0, 0), &grid.W(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.P(0, 0, 0), &grid.P(-1, j, k), sizeof(Real) * grid.structure.gx);
    };

    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        MPI_Isend(bufferOut.data, cast_int(bufferOut.node_dim) * 4, Real_MPI,
                  neigh_rank, SOUTH_EAST_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        MPI_Irecv(bufferIn.data, cast_int(bufferIn.node_dim) * 4, Real_MPI,
                  neigh_rank, NORTH_WEST_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);
        const index_t j = -1;
        const index_t k = grid.structure.nz;
        memcpy(&grid.U(-1, j, k), &buffer.U(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.V(-1, j, k), &buffer.V(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.W(-1, j, k), &buffer.W(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.P(-1, j, k), &buffer.P(0, 0, 0), sizeof(Real) * grid.structure.gx);
    };
}

/// SOUTH WEST /////////////////////////////////////////////////////////////////////////////////////////
namespace shwt {
    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t j = 0;
        const index_t k = 0;
        memcpy(&bufferOut.U(0, 0, 0), &grid.U(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.V(0, 0, 0), &grid.V(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.W(0, 0, 0), &grid.W(-1, j, k), sizeof(Real) * grid.structure.gx);
        memcpy(&bufferOut.P(0, 0, 0), &grid.P(-1, j, k), sizeof(Real) * grid.structure.gx);
    };

    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        MPI_Isend(bufferOut.data, cast_int(bufferOut.node_dim) * 4, Real_MPI,
                  neigh_rank, SOUTH_WEST_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        MPI_Irecv(bufferIn.data, cast_int(bufferIn.node_dim) * 4, Real_MPI,
                  neigh_rank, NORTH_EAST_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);
        const index_t j = -1;
        const index_t k = -1;
        memcpy(&grid.U(-1, j, k), &buffer.U(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.V(-1, j, k), &buffer.V(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.W(-1, j, k), &buffer.W(0, 0, 0), sizeof(Real) * grid.structure.gx);
        memcpy(&grid.P(-1, j, k), &buffer.P(0, 0, 0), sizeof(Real) * grid.structure.gx);
    };
}


/**
 * Builds boundary conditions for a unique space (unused)
 */
inline void buildBoundaries(Boundaries &boundaries, const boundaryDomainFunctions &boundaryFunctions) {
    Condition *northCond, *southCond, *eastCond, *westCond, *frontCond, *backCond;

    /// Define face mappers ////////////////////////////////////////////////////////////////////////////////
    // North
    PhysicalCondition::Mapper northFace = [](GridData &grid,
                                             const Real currentTime,
                                             const boundaryFaceFunctions &functions) {
        getStaggeredSpacing(grid, sdx, sdy, sdz);
        getExactFunctions(functions, eU, eV, eW);

        const index_t j = grid.structure.ny;
        const Real y = real(j + grid.structure.py) * grid.structure.dy;

        // apply on face with y constant

        for (index_t k = 0; k < grid.structure.nz; k++)
            for (index_t i = 0; i < grid.structure.nx; i++) {
                Real x = real(i + grid.structure.px) * grid.structure.dx;
                Real z = real(k + grid.structure.pz) * grid.structure.dz;

                // On y = phy_dim for domain point we have exact for V
                grid.V(i, j - 1, k) = eV(x + sdx, y, z + sdz, currentTime);
                // For ghost points we have useless V, other approximate

                grid.U(i, j, k) = 2 * eU(x + grid.structure.dx, y, z + sdz, currentTime)
                                  - grid.U(i, j - 1, k);
                grid.V(i, j, k) = 0;
                grid.W(i, j, k) = 2 * eW(x + sdx, y, z + grid.structure.dz, currentTime)
                                  - grid.W(i, j - 1, k);
            }
    };

    northCond = new PhysicalCondition(northFace, boundaryFunctions[NORTH_FACE_ID]);
    boundaries.addCond(*northCond);

    // South
    PhysicalCondition::Mapper southFace = [](GridData &grid,
                                             const Real currentTime,
                                             const boundaryFaceFunctions &functions) {
        getStaggeredSpacing(grid, sdx, sdy, sdz);
        getExactFunctions(functions, eU, eV, eW);

        const index_t j = 0;
        const Real y = real(j + grid.structure.py) * grid.structure.dy;

        // apply on face with y constant

        for (index_t k = 0; k < grid.structure.nz; k++)
            for (index_t i = 0; i < grid.structure.nx; i++) {
                Real x = real(i + grid.structure.px) * grid.structure.dx;
                Real z = real(k + grid.structure.pz) * grid.structure.dz;

                // On y = 0 for ghost point we hae exact for V, other approximate
                grid.U(i, j - 1, k) = 2 * eU(x + grid.structure.dx, y, z + sdz, currentTime)
                                      - grid.U(i, j, k);
                grid.V(i, j - 1, k) = eV(x + sdx, y, z + sdz, currentTime);
                grid.W(i, j - 1, k) = 2 * eW(x + sdx, y, z + grid.structure.dz, currentTime)
                                      - grid.W(i, j, k);
            }
    };

    southCond = new PhysicalCondition(southFace, boundaryFunctions[SOUTH_FACE_ID]);
    boundaries.addCond(*southCond);

    // East
    PhysicalCondition::Mapper eastFace = [](GridData &grid,
                                            const Real currentTime,
                                            const boundaryFaceFunctions &functions) {
        getStaggeredSpacing(grid, sdx, sdy, sdz);
        getExactFunctions(functions, eU, eV, eW);

        const index_t k = grid.structure.nz;
        const Real z = real(k + grid.structure.pz) * grid.structure.dz;

        for (index_t j = 0; j < grid.structure.ny; j++)
            for (index_t i = 0; i < grid.structure.nx; i++) {
                Real x = real(i + grid.structure.px) * grid.structure.dx;
                Real y = real(j + grid.structure.py) * grid.structure.dy;

                // On z = phy_dim for domain point we have exact for W
                grid.W(i, j, k - 1) = eW(x + sdx, y + sdy, z, currentTime);
                // For ghost points we gave useless W, other interpolate
                grid.U(i, j, k) = 2 * eU(x + grid.structure.dx, y + sdy, z, currentTime)
                                  - grid.U(i, j, k - 1);
                grid.V(i, j, k) = 2 * eV(x + sdx, y + grid.structure.dy, z, currentTime)
                                  - grid.V(i, j, k - 1);
                grid.W(i, j, k) = 0;
            }
    };

    eastCond = new PhysicalCondition(eastFace, boundaryFunctions[EAST_FACE_ID]);
    boundaries.addCond(*eastCond);

    // West
    PhysicalCondition::Mapper westFace = [](GridData &grid,
                                            const Real currentTime,
                                            const boundaryFaceFunctions &functions) {
        getStaggeredSpacing(grid, sdx, sdy, sdz);
        getExactFunctions(functions, eU, eV, eW);

        const index_t k = 0;
        const Real z = real(k + grid.structure.pz) * grid.structure.dz;

        for (index_t j = 0; j < grid.structure.ny; j++)
            for (index_t i = 0; i < grid.structure.nx; i++) {
                Real x = real(i + grid.structure.px) * grid.structure.dx;
                Real y = real(j + grid.structure.py) * grid.structure.dy;

                // On z = 0 for ghost point we have exact for W, other approximate
                grid.U(i, j, k - 1) = 2 * eU(x + grid.structure.dx, y + sdy, z, currentTime)
                                      - grid.U(i, j, k);
                grid.V(i, j, k - 1) = 2 * eV(x + sdx, y + grid.structure.dy, z, currentTime)
                                      - grid.V(i, j, k);
                grid.W(i, j, k - 1) = eW(x + sdx, y + sdy, z, currentTime);
            }
    };

    westCond = new PhysicalCondition(westFace, boundaryFunctions[WEST_FACE_ID]);
    boundaries.addCond(*westCond);

    // Front
    PhysicalCondition::Mapper frontFace = [](GridData &grid,
                                             const Real currentTime,
                                             const boundaryFaceFunctions &functions) {
        getStaggeredSpacing(grid, sdx, sdy, sdz);
        getExactFunctions(functions, eU, eV, eW);

        const index_t i = 0;
        const Real x = real(i + grid.structure.px) * grid.structure.dx;

        for (index_t k = 0; k < grid.structure.nz; k++)
            for (index_t j = 0; j < grid.structure.ny; j++) {
                Real y = real(j + grid.structure.py) * grid.structure.dy;
                Real z = real(k + grid.structure.pz) * grid.structure.dz;

                // On x = 0 for ghost point we have exact for U, other approximate
                grid.U(i - 1, j, k) = eU(x, y + sdy, z + sdz, currentTime);
                grid.V(i - 1, j, k) = 2 * eV(x, y + grid.structure.dy, z + sdz, currentTime)
                                      - grid.V(i, j, k);
                grid.W(i - 1, j, k) = 2 * eW(x, y + sdy, z + grid.structure.dz, currentTime)
                                      - grid.W(i, j, k);
            }
    };

    frontCond = new PhysicalCondition(frontFace, boundaryFunctions[FRONT_FACE_ID]);
    boundaries.addCond(*frontCond);

    // Back
    PhysicalCondition::Mapper backFace = [](GridData &grid,
                                            const Real currentTime,
                                            const boundaryFaceFunctions &functions) {
        getStaggeredSpacing(grid, sdx, sdy, sdz);
        getExactFunctions(functions, eU, eV, eW);

        const index_t i = grid.structure.nx;
        const Real x = real(i + grid.structure.px) * grid.structure.dx;

        for (index_t k = 0; k < grid.structure.nz; k++)
            for (index_t j = 0; j < grid.structure.ny; j++) {
                Real y = real(j + grid.structure.py) * grid.structure.dy;
                Real z = real(k + grid.structure.pz) * grid.structure.dz;

                // On x = phy_dim for domain point we have exact for U
                grid.U(i - 1, j, k) = eU(x, y + sdy, z + sdz, currentTime);
                // For ghost point we have useless U, other approximate
                grid.U(i, j, k) = 0;
                grid.V(i, j, k) = 2 * eV(x, y + grid.structure.dy, z + sdz, currentTime)
                                  - grid.V(i - 1, j, k);
                grid.W(i, j, k) = 2 * eW(x, y + sdy, z + grid.structure.dz, currentTime)
                                  - grid.W(i - 1, j, k);
            }
    };

    backCond = new PhysicalCondition(backFace, boundaryFunctions[BACK_FACE_ID]);
    boundaries.addCond(*backCond);
}

/**
 * Builds boundary conditions for a partitioned space
 */
inline void buildMPIBoundaries(const C2Decomp &decomp, const GridStructure &gridStructure, MPIBoundaries &boundaries,
                               const boundaryDomainFunctions &boundaryFunctions) {
    /// Determine where the domain is positioned ///////////////////////////////////////////////////////////
    // global position of this process
    int n_y_proc = decomp.dims[0];
    int n_z_proc = decomp.dims[1];
    int this_y_pos = decomp.coord[0];
    int this_z_pos = decomp.coord[1];

    // flags to define if the processor is on a physical boundary
    bool isOnTop = (this_y_pos == n_y_proc - 1);
    bool isOnBottom = (this_y_pos == 0);
    bool isOnLeft = (this_z_pos == 0);
    bool isOnRight = (this_z_pos == n_z_proc - 1);

    /// Define face mappers ////////////////////////////////////////////////////////////////////////////////

    // NORTH (MAX Y)
    if (isOnTop) {
        // The processor is on top of global domain, it has to apply on the upper face the physical BC
        // Build the #condition that will associate the mapper (#face) to the corresponding #BC_functions
        auto *northCond = new PhysicalCondition(north::face, boundaryFunctions[NORTH_FACE_ID]);
        // Add the #condition to the collection
        boundaries.addCond(*northCond);
    } else {
        // The processor is not on top of the global domain,
        // so we have to set up a communication layer with the above processor

        // process rank that is to the north of this one
        const int north_neigh_rank = decomp.neighbor[0][2];

        // Define the structure of communication buffers,
        // Since north face is a X-Z plane the face has dimension 1 for Y
        auto *bufferStructure = new GridStructure({gridStructure.gx, 1, gridStructure.gz}, {0, 0, 0}, {0, 0, 0}, 0);

        // Create the BC, I assign as neighbour the north_neighbour
        auto *northCond = new MPICondition(north::init, north::exc, north::mapp, *bufferStructure, north_neigh_rank);

        // Add the BC to the collection
        boundaries.addMPICond(*northCond);
    }

    // SOUTH (MIN Y)
    if (isOnBottom) {
        // The processor is on bottom of global domain, it has to apply on the lower face the physical BC
        auto southCond = new PhysicalCondition(south::face, boundaryFunctions[SOUTH_FACE_ID]);
        boundaries.addCond(*southCond);
    } else {
        // The processor is not on bottom of the global domain,
        // so we have to set up a communication layer with the below processor

        // process rank that is to the south
        const int south_neigh_rank = decomp.neighbor[0][3];

        // Since south face is a X-Z plane the face has dimension 1 for Y
        auto *bufferStructure = new GridStructure({gridStructure.gx, 1, gridStructure.gz}, {0, 0, 0}, {0, 0, 0}, 0);
        auto *southCond = new MPICondition(south::init, south::exc, south::mapp, *bufferStructure, south_neigh_rank);
        boundaries.addMPICond(*southCond);
    }

    // EAST (MAX Z)
    if (isOnRight) {
        // The processor is at the right of global domain, it has to apply on the right face the physical BC
        auto eastCond = new PhysicalCondition(east::face, boundaryFunctions[EAST_FACE_ID]);
        boundaries.addCond(*eastCond);
    } else {
        // The processor is not on the right of the global domain,
        // so we have to set up a communication layer with the processor at right

        // process rank that is to the east
        const int east_neigh_rank = decomp.neighbor[0][4];

        // Since east face is a X-Y plane the face has dimension 1 for Z
        auto *bufferStructure = new GridStructure({gridStructure.gx, gridStructure.gy, 1}, {0, 0, 0}, {0, 0, 0}, 0);
        auto *eastCond = new MPICondition(east::init, east::exc, east::mapp, *bufferStructure, east_neigh_rank);
        boundaries.addMPICond(*eastCond);
    }

    // WEST (MIN Z)
    if (isOnLeft) {
        // The processor is on left of global domain, it has to apply on the left face the physical BC
        auto westCond = new PhysicalCondition(west::face, boundaryFunctions[WEST_FACE_ID]);
        boundaries.addCond(*westCond);
    } else {
        // The processor is not on the left of the global domain,
        // so we have to set up a communication layer with the processor at left

        // process rank that is to the west
        const int west_proc_rank = decomp.neighbor[0][5];

        // Since east face is a X-Y plane the face has dimension 1 for Z
        auto *bufferStructure = new GridStructure({gridStructure.gx, gridStructure.gy, 1}, {0, 0, 0}, {0, 0, 0}, 0);
        auto *westCond = new MPICondition(west::init, west::exc, west::mapp, *bufferStructure, west_proc_rank);
        boundaries.addMPICond(*westCond);
    }


    if (!(isOnTop || isOnRight)) {
        int nhet_proc_rank;
        const int coord[] = {this_y_pos + 1, this_z_pos + 1};
        MPI_Cart_rank(decomp.DECOMP_2D_COMM_CART_X, coord, &nhet_proc_rank);

        auto *bufferStructure = new GridStructure({gridStructure.gx, 1, 1}, {0, 0, 0}, {0, 0, 0}, 0);
        auto *nhetCond = new MPICondition(nhet::init, nhet::exc, nhet::mapp, *bufferStructure, nhet_proc_rank);
        boundaries.addMPICond(*nhetCond);
    }

    if (!(isOnTop || isOnLeft)) {
        int nhwt_proc_rank;
        const int coord[] = {this_y_pos + 1, this_z_pos - 1};
        MPI_Cart_rank(decomp.DECOMP_2D_COMM_CART_X, coord, &nhwt_proc_rank);


        auto *bufferStructure = new GridStructure({gridStructure.gx, 1, 1}, {0, 0, 0}, {0, 0, 0}, 0);
        auto *nhwtCond = new MPICondition(nhwt::init, nhwt::exc, nhwt::mapp, *bufferStructure, nhwt_proc_rank);
        boundaries.addMPICond(*nhwtCond);
    }

    if (!(isOnBottom || isOnRight)) {
        int shet_proc_rank;
        const int coord[] = {this_y_pos - 1, this_z_pos + 1};
        MPI_Cart_rank(decomp.DECOMP_2D_COMM_CART_X, coord, &shet_proc_rank);

        auto *bufferStructure = new GridStructure({gridStructure.gx, 1, 1}, {0, 0, 0}, {0, 0, 0}, 0);
        auto *shetCond = new MPICondition(shet::init, shet::exc, shet::mapp, *bufferStructure, shet_proc_rank);
        boundaries.addMPICond(*shetCond);
    }

    if (!(isOnBottom || isOnLeft)) {
        int shwt_proc_rank;
        const int coord[] = {this_y_pos - 1, this_z_pos - 1};
        MPI_Cart_rank(decomp.DECOMP_2D_COMM_CART_X, coord, &shwt_proc_rank);

        auto *bufferStructure = new GridStructure({gridStructure.gx, 1, 1}, {0, 0, 0}, {0, 0, 0}, 0);
        auto *shwtCond = new MPICondition(shwt::init, shwt::exc, shwt::mapp, *bufferStructure, shwt_proc_rank);
        boundaries.addMPICond(*shwtCond);
    }

    // Font and back faces are always physical boundaries in pencil domain decomposition

    // FRONT (MIN X)
    auto frontCond = new PhysicalCondition(front::face, boundaryFunctions[FRONT_FACE_ID]);
    boundaries.addCond(*frontCond);

    // BACK (MAX X)
    auto backCond = new PhysicalCondition(back::face, boundaryFunctions[BACK_FACE_ID]);
    boundaries.addCond(*backCond);
}
