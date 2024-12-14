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

/// NORTH //////////////////////////////////////////////////////////////////////////////////////////////
namespace north {
    // This lambda defines how the #functions have to be applied to the #grid
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const std::vector<TFunction> &functions) {

        // Use macro to get some variables
        getStaggeredSpacing(grid, sdx, sdy, sdz);
        getExactFunctions(functions, eU, eV, eW);

        // The upper boundary is at the max Y
        const index_t j = grid.structure.ny;
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
            }
    };
    // This lambda defines how the outgoing communication buffer has to be initialized
    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t j = grid.structure.ny - 1;
        for (index_t k = 0; k < grid.structure.nz; k++) {
            memcpy(&bufferOut.U(0, 0, k), &grid.U(0, j, k), sizeof(Real) * grid.structure.nx);
            memcpy(&bufferOut.V(0, 0, k), &grid.V(0, j, k), sizeof(Real) * grid.structure.nx);
            memcpy(&bufferOut.W(0, 0, k), &grid.W(0, j, k), sizeof(Real) * grid.structure.nx);
        }
    };

    // This lambda defines which data has to be shared with the given neighbour (#neigh_rank)
    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        // This proc will send his outgoing buffer with the tag #NORTH_BUFFER_TAG
        // (means that the buffer is the top layer of this domain)
        MPI_Isend(bufferOut.velocity_data, int(bufferOut.node_dim) * 3, Real_MPI,
                  neigh_rank, NORTH_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        // This proc will receive into ingoing buffer data with tag #SOUTH_BUFFER_TAG
        // (means that the buffer is the top layer of neighbour domain)
        MPI_Irecv(bufferIn.velocity_data, int(bufferIn.node_dim) * 3, Real_MPI,
                  neigh_rank, SOUTH_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    // This lambda defines how the ingoing communication buffer has to be copied into the local model grid
    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        // Since communication are asynchronous I wait for the ingoing buffer to be completely received
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);

        // Copy the buffer into the upper ghost point layer
        const index_t j = grid.structure.ny;
        for (index_t k = 0; k < grid.structure.nz; k++) {
            memcpy(&grid.U(0, j, k), &buffer.U(0, 0, k), sizeof(Real) * grid.structure.nx);
            memcpy(&grid.V(0, j, k), &buffer.V(0, 0, k), sizeof(Real) * grid.structure.nx);
            memcpy(&grid.W(0, j, k), &buffer.W(0, 0, k), sizeof(Real) * grid.structure.nx);
        }
    };
}
/// SOUTH //////////////////////////////////////////////////////////////////////////////////////////////
namespace south {
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const std::vector<TFunction> &functions) {

        getStaggeredSpacing(grid, sdx, sdy, sdz);
        getExactFunctions(functions, eU, eV, eW);

        // South face has y=0
        const index_t j = 0;
        const Real y = real(j + grid.structure.py) * grid.structure.dy;

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

    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t j = 0;
        for (index_t k = 0; k < grid.structure.nz; k++) {
            memcpy(&bufferOut.U(0, 0, k), &grid.U(0, j, k), sizeof(Real) * grid.structure.nx);
            memcpy(&bufferOut.V(0, 0, k), &grid.V(0, j, k), sizeof(Real) * grid.structure.nx);
            memcpy(&bufferOut.W(0, 0, k), &grid.W(0, j, k), sizeof(Real) * grid.structure.nx);
        }
    };

    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        MPI_Isend(bufferOut.velocity_data, int(bufferOut.node_dim) * 3, Real_MPI,
                  neigh_rank, SOUTH_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        MPI_Irecv(bufferIn.velocity_data, int(bufferIn.node_dim) * 3, Real_MPI,
                  neigh_rank, NORTH_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);

        // Copy the buffer into the ghost point layer
        const index_t j = -1;
        for (index_t k = 0; k < grid.structure.nz; k++) {
            memcpy(&grid.U(0, j, k), &buffer.U(0, 0, k), sizeof(Real) * grid.structure.nx);
            memcpy(&grid.V(0, j, k), &buffer.V(0, 0, k), sizeof(Real) * grid.structure.nx);
            memcpy(&grid.W(0, j, k), &buffer.W(0, 0, k), sizeof(Real) * grid.structure.nx);
        }
    };
}
/// EAST ///////////////////////////////////////////////////////////////////////////////////////////////
namespace east {
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const std::vector<TFunction> &functions) {

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

    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t k = grid.structure.nz - 1;

        for (index_t j = 0; j < grid.structure.ny; j++) {
            memcpy(&bufferOut.U(0, j, 0), &grid.U(0, j, k), sizeof(Real) * grid.structure.nx);
            memcpy(&bufferOut.V(0, j, 0), &grid.V(0, j, k), sizeof(Real) * grid.structure.nx);
            memcpy(&bufferOut.W(0, j, 0), &grid.W(0, j, k), sizeof(Real) * grid.structure.nx);
        }
    };

    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        MPI_Isend(bufferOut.velocity_data, int(bufferOut.node_dim) * 3, Real_MPI, neigh_rank,
                  EAST_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        MPI_Irecv(bufferIn.velocity_data, int(bufferIn.node_dim) * 3, Real_MPI, neigh_rank,
                  WEST_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);

        // Copy the buffer into the ghost point layer
        const index_t k = grid.structure.nz;
        for (index_t j = 0; j < grid.structure.ny; j++) {
            memcpy(&grid.U(0, j, k), &buffer.U(0, j, 0), sizeof(Real) * grid.structure.nx);
            memcpy(&grid.V(0, j, k), &buffer.V(0, j, 0), sizeof(Real) * grid.structure.nx);
            memcpy(&grid.W(0, j, k), &buffer.W(0, j, 0), sizeof(Real) * grid.structure.nx);
        }
    };
}
/// WEST ///////////////////////////////////////////////////////////////////////////////////////////////
namespace west {
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const std::vector<TFunction> &functions) {

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

    MPICondition::BufferInitializer init = [](GridData &grid, GridData &bufferOut) {
        // I want to copy the last in-domain layer
        const index_t k = 0;
        for (index_t j = 0; j < grid.structure.ny; j++) {
            memcpy(&bufferOut.U(0, j, 0), &grid.U(0, j, k), sizeof(Real) * grid.structure.nx);
            memcpy(&bufferOut.V(0, j, 0), &grid.V(0, j, k), sizeof(Real) * grid.structure.nx);
            memcpy(&bufferOut.W(0, j, 0), &grid.W(0, j, k), sizeof(Real) * grid.structure.nx);
        }
    };

    MPICondition::BufferExchanger exc = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut,
                                           MPI_Request *requestIn, int neigh_rank) {
        MPI_Isend(bufferOut.velocity_data, int(bufferOut.node_dim) * 3, Real_MPI, neigh_rank,
                  WEST_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        MPI_Irecv(bufferIn.velocity_data, int(bufferIn.node_dim) * 3, Real_MPI, neigh_rank,
                  EAST_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    MPICondition::BufferMapper mapp = [](GridData &grid, GridData &buffer, MPI_Request *requestOut,
                                         MPI_Request *requestIn) {
        MPI_Wait(requestIn, MPI_STATUS_IGNORE);

        // Copy the buffer into the ghost point layer
        const index_t k = -1;
        for (index_t j = 0; j < grid.structure.ny; j++) {
            memcpy(&grid.U(0, j, k), &buffer.U(0, j, 0), sizeof(Real) * grid.structure.nx);
            memcpy(&grid.V(0, j, k), &buffer.V(0, j, 0), sizeof(Real) * grid.structure.nx);
            memcpy(&grid.W(0, j, k), &buffer.W(0, j, 0), sizeof(Real) * grid.structure.nx);
        }
    };
}
/// FRONT //////////////////////////////////////////////////////////////////////////////////////////////
namespace front {
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const std::vector<TFunction> &functions) {

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
}
/// BACK ///////////////////////////////////////////////////////////////////////////////////////////////
namespace back {
    PhysicalCondition::Mapper face = [](GridData &grid,
                                        const Real currentTime,
                                        const std::vector<TFunction> &functions) {

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
}


/**
 * Builds boundary conditions for a unique space (unused)
 */
inline void buildBoundaries(Boundaries &boundaries, const std::vector<TFunction> &boundaryFunctions) {
    Condition *northCond, *southCond, *eastCond, *westCond, *frontCond, *backCond;

    /// Define face mappers ////////////////////////////////////////////////////////////////////////////////
    // North
    PhysicalCondition::Mapper northFace = [](GridData &grid,
                                             const Real currentTime,
                                             const std::vector<TFunction> &functions) {

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

    northCond = new PhysicalCondition(northFace, boundaryFunctions);
    boundaries.addCond(*northCond);

    // South
    PhysicalCondition::Mapper southFace = [](GridData &grid,
                                             const Real currentTime,
                                             const std::vector<TFunction> &functions) {

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

    southCond = new PhysicalCondition(southFace, boundaryFunctions);
    boundaries.addCond(*southCond);

    // East
    PhysicalCondition::Mapper eastFace = [](GridData &grid,
                                            const Real currentTime,
                                            const std::vector<TFunction> &functions) {

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

    eastCond = new PhysicalCondition(eastFace, boundaryFunctions);
    boundaries.addCond(*eastCond);

    // West
    PhysicalCondition::Mapper westFace = [](GridData &grid,
                                            const Real currentTime,
                                            const std::vector<TFunction> &functions) {

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

    westCond = new PhysicalCondition(westFace, boundaryFunctions);
    boundaries.addCond(*westCond);

    // Front
    PhysicalCondition::Mapper frontFace = [](GridData &grid,
                                             const Real currentTime,
                                             const std::vector<TFunction> &functions) {

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

    frontCond = new PhysicalCondition(frontFace, boundaryFunctions);
    boundaries.addCond(*frontCond);

    // Back
    PhysicalCondition::Mapper backFace = [](GridData &grid,
                                            const Real currentTime,
                                            const std::vector<TFunction> &functions) {

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

    backCond = new PhysicalCondition(backFace, boundaryFunctions);
    boundaries.addCond(*backCond);
}

/**
 * Builds boundary conditions for a partitioned space
 */
inline void buildMPIBoundaries(const C2Decomp &decomp, const GridStructure &gridStructure, MPIBoundaries &boundaries,
                               const std::vector<TFunction> &boundaryFunctions) {

#define int(v) static_cast<int>(v)

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
        auto *northCond = new PhysicalCondition(north::face, boundaryFunctions);
        // Add the #condition to the collection
        boundaries.addCond(*northCond);
    } else {
        // The processor is not on top of the global domain,
        // so we have to set up a communication layer with the above processor

        // process rank that is to the north of this one
        const int north_neigh_rank = decomp.neighbor[0][2];

        // Define the structure of communication buffers,
        // Since north face is a X-Z plane the face has dimension 1 for Y
        auto *bufferStructure = new GridStructure({gridStructure.nx, 1, gridStructure.nz}, {0, 0, 0}, {0, 0, 0}, 0);

        // Create the BC, I assign as neighbour the north_neighbour
        auto *northCond = new MPICondition(north::init, north::exc, north::mapp, *bufferStructure, north_neigh_rank);

        // Add the BC to the collection
        boundaries.addMPICond(*northCond);
    }

    // SOUTH (MIN Y)
    if (isOnBottom) {
        // The processor is on bottom of global domain, it has to apply on the lower face the physical BC
        auto southCond = new PhysicalCondition(south::face, boundaryFunctions);
        boundaries.addCond(*southCond);
    } else {
        // The processor is not on bottom of the global domain,
        // so we have to set up a communication layer with the below processor

        // process rank that is to the south
        const int south_neigh_rank = decomp.neighbor[0][3];

        // Since south face is a X-Z plane the face has dimension 1 for Y
        auto *bufferStructure = new GridStructure({gridStructure.nx, 1, gridStructure.nz}, {0, 0, 0}, {0, 0, 0}, 0);
        auto *southCond = new MPICondition(south::init, south::exc, south::mapp, *bufferStructure, south_neigh_rank);
        boundaries.addMPICond(*southCond);
    }

    // EAST (MAX Z)
    if (isOnRight) {
        // The processor is at the right of global domain, it has to apply on the right face the physical BC
        auto eastCond = new PhysicalCondition(east::face, boundaryFunctions);
        boundaries.addCond(*eastCond);
    } else {
        // The processor is not on the right of the global domain,
        // so we have to set up a communication layer with the processor at right

        // process rank that is to the east
        const int east_neigh_rank = decomp.neighbor[0][4];

        // Since east face is a X-Y plane the face has dimension 1 for Z
        auto *bufferStructure = new GridStructure({gridStructure.nx, gridStructure.ny, 1}, {0, 0, 0}, {0, 0, 0}, 0);
        auto *eastCond = new MPICondition(east::init, east::exc, east::mapp, *bufferStructure, east_neigh_rank);
        boundaries.addMPICond(*eastCond);
    }

    // WEST (MIN Z)
    if (isOnLeft) {
        // The processor is on left of global domain, it has to apply on the left face the physical BC
        auto westCond = new PhysicalCondition(west::face, boundaryFunctions);
        boundaries.addCond(*westCond);
    } else {
        // The processor is not on the left of the global domain,
        // so we have to set up a communication layer with the processor at left

        // process rank that is to the west
        const int west_proc_rank = decomp.neighbor[0][5];

        // Since east face is a X-Y plane the face has dimension 1 for Z
        auto *bufferStructure = new GridStructure({gridStructure.nx, gridStructure.ny, 1}, {0, 0, 0}, {0, 0, 0}, 0);
        auto *westCond = new MPICondition(west::init, west::exc, west::mapp, *bufferStructure, west_proc_rank);
        boundaries.addMPICond(*westCond);
    }



    // Font and back faces are always physical boundaries in pencil domain decomposition

    // FRONT (MIN X)
    auto frontCond = new PhysicalCondition(front::face, boundaryFunctions);
    boundaries.addCond(*frontCond);

    // BACK (MAX X)
    auto backCond = new PhysicalCondition(back::face, boundaryFunctions);
    boundaries.addCond(*backCond);
}
