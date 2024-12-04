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

/**
 * Builds boundary conditions for a unique space
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


inline void buildMPIBoundaries(const C2Decomp &decomp, const GridStructure &gridStructure, MPIBoundaries &boundaries,
                               const std::vector<TFunction> &boundaryFunctions) {

    /// Determine where the domain is positioned ///////////////////////////////////////////////////////////
    // global reference of this process position
    int n_y_proc = decomp.dims[0] - 1;
    int n_z_proc = decomp.dims[1] - 1;
    int this_y_pos = decomp.coord[0];
    int this_z_pos = decomp.coord[1];

    // flags to define if the processor is on a physical boundary
    bool isOnTop = (this_y_pos == n_y_proc - 1);
    bool isOnBottom = (this_y_pos == 0);
    bool isOnLeft = (this_z_pos == 0);
    bool isOnRight = (this_z_pos == n_z_proc - 1);

    /// Define face mappers ////////////////////////////////////////////////////////////////////////////////

    // the buffer exchanger will be the same for all
    MPICondition::BufferExchanger exchanger = [](GridData &bufferOut, GridData &bufferIn, MPI_Request *requestIn, MPI_Request *requestOut,
                                                 int proc_rank) {
        MPI_Isend(bufferOut.velocity_data, bufferOut.dim * 3, Real_MPI, proc_rank, SOUTH_BUFFER_TAG, MPI_COMM_WORLD, requestOut);
        MPI_Irecv(bufferIn.velocity_data, bufferIn.dim * 3, Real_MPI, proc_rank, NORTH_BUFFER_TAG, MPI_COMM_WORLD, requestIn);
    };

    // North
    if (isOnTop) {
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

        auto *northCond = new PhysicalCondition(northFace, boundaryFunctions);
        boundaries.addCond(*northCond);

    } else {
        // TODO DETERMINE TO WHICH PROCESSOR I HAVE TO COMMUNICATE

        // process rank that is to the north
        const int proc_rank = decomp.neighbor[0][5];

        MPICondition::BufferInitializer northInit = [](GridData &grid, GridData &bufferOut) {
            // I want to copy the last in-domain layer
            const index_t j = grid.structure.ny - 1;

            // Be careful when working with buffer since its always on x-z plane
            for (index_t k = 0; k < grid.structure.nz; k++) {
                memcpy(&bufferOut.U(0, 0, k), &grid.U(0, j, k), sizeof(Real) * grid.structure.nx);
                memcpy(&bufferOut.V(0, 0, k), &grid.V(0, j, k), sizeof(Real) * grid.structure.nx);
                memcpy(&bufferOut.W(0, 0, k), &grid.W(0, j, k), sizeof(Real) * grid.structure.nx);
            }
        };

        MPICondition::BufferMapper northMapp = [](GridData &grid, GridData &buffer, MPI_Request *requestIn, MPI_Request *requestOut) {
            MPI_Status status{};
            MPI_Wait(requestIn, &status);

            MPI_Request_free(requestIn);
            MPI_Request_free(requestOut);

            // Copy the buffer into the ghost point layer
            const index_t j = grid.structure.ny;
            for (index_t k = 0; k < grid.structure.nz; k++) {
                memcpy(&grid.U(0, j, k), &buffer.U(0, 0, k), sizeof(Real) * grid.structure.nx);
                memcpy(&grid.V(0, j, k), &buffer.V(0, 0, k), sizeof(Real) * grid.structure.nx);
                memcpy(&grid.W(0, j, k), &buffer.W(0, 0, k), sizeof(Real) * grid.structure.nx);
            }
        };

        // Since north face is a X-Z plane the face has dimension 1 for Y
        GridStructure bufferStructure({gridStructure.nx, 1, gridStructure.nz}, {0,0,0}, {0,0,0}, 0);
        auto *northCond = new MPICondition(northInit, exchanger, northMapp, bufferStructure, proc_rank);

        // NOTE that we add mpi condition for both precondition and condition
        boundaries.addPreCond(*northCond).addCond(*northCond);
    }

    //TODO SOUTH
    if (isOnBottom) {

    } else {
        // process rank that is to the south
        const int proc_rank = decomp.neighbor[0][4];
    }

    //TODO EAST
    if (isOnRight) {

    } else {
        // process rank that is to the east
        const int proc_rank = decomp.neighbor[0][3];
    }

    //TODO WEST
    if (isOnLeft) {

    } else {
        // process rank that is to the west
        const int proc_rank = decomp.neighbor[0][2];
    }

    //TODO FRONT

    //TODO BACK


}