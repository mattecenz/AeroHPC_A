#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "Traits.hpp"
#include "C2Decomp.hpp"

class Parameters {
public:
    // Spatial Info
    Real dimX, dimY, dimZ; // Physical dimensions of global space
    Real originX, originY, originZ; // Physical origin of global space

    // Logical Grid Info
    index_t glob_nX, glob_nY, glob_nZ; // Number of global nodes
    index_t loc_nX, loc_nY, loc_nZ; // Number of local nodes
    index_t st_nX, st_nY, st_nZ; // Local nodes starting global index
    index_t loc_gnX, loc_gnY, loc_gnZ; // Number of local nodes + ghosts
    index_t grid_ndim; // Total number of local nodes
    index_t grid_gndim; // Total number of local nodes + ghosts

    // Physical Grid Info
    Real dX, dY, dZ; // Physical spacing between nodes
    Real dX2, dY2, dZ2; // Physical spacing between staggered nodes
    index_t phy_nX, phy_nY, phy_nZ; // Physical grid number of nodes (1 more than logical)
    index_t phy_ndim; // Physical grid total number of nodes

    // SubDomain domain info (rank of neighbours)
    int neigh_north, neigh_south,
            neigh_east, neigh_west,
            neigh_north_east, neigh_north_west,
            neigh_south_east, neigh_south_west;

    // SubDomain boundary Info
    bool isOnTop, isOnBottom, isOnRight, isOnLeft; // Where is globally positioned
    bool periodicX, periodicY, periodicZ; // If the boundary is periodic on directions
    bool hasPeriodicLayerX, hasPeriodicLayerY, hasPeriodicLayerZ;

    // MPI Types
    MPI_Datatype XYFace;
    MPI_Datatype XZFace;
    MPI_Datatype XRow;

    MPI_Datatype XYFace_Periodic;
    MPI_Datatype XZFace_Periodic;

    // Time Info
    Real dt; // Physical delta Time between timesteps

    // Solver Info
    index_t timesteps; // Number of timesteps to perform
    Real Re; // Reynolds number

    Parameters(const Real dimX, const Real dimY, const Real dimZ,
               const Real originX, const Real originY, const Real originZ,
               const Real dt, const index_t timesteps, const Real Re,
               const bool periodicBC[3], const C2Decomp &c2D)
        : dimX(dimX), dimY(dimY), dimZ(dimZ),
          originX(originX), originY(originY), originZ(originZ),
          dt(dt), timesteps(timesteps), Re(Re) {

        periodicX = periodicBC[0];
        periodicY = periodicBC[1];
        periodicZ = periodicBC[2];

        // global position of this process
        const int n_y_proc = c2D.dims[0];
        const int n_z_proc = c2D.dims[1];
        const int this_y_pos = c2D.coord[0];
        const int this_z_pos = c2D.coord[1];

        // flags to define if the processor is on a physical boundary
        isOnTop = (this_y_pos == n_y_proc - 1);
        isOnBottom = (this_y_pos == 0);
        isOnLeft = (this_z_pos == 0);
        isOnRight = (this_z_pos == n_z_proc - 1);

        hasPeriodicLayerX = periodicX;
        hasPeriodicLayerY = periodicY && isOnBottom;
        hasPeriodicLayerZ = periodicZ && isOnRight;

        glob_nX = c2D.nxGlobal;
        glob_nY = c2D.nyGlobal;
        glob_nZ = c2D.nzGlobal;

        dX = dimX / real(glob_nX - periodicX);
        dY = dimY / real(glob_nY - periodicY);
        dZ = dimZ / real(glob_nZ - periodicZ);

        dX2 = dX / real(2);
        dY2 = dY / real(2);
        dZ2 = dZ / real(2);

        loc_nX = c2D.xSize[0];
        loc_nY = c2D.xSize[1];
        loc_nZ = c2D.xSize[2];

        st_nX = c2D.xStart[0];
        st_nY = c2D.xStart[1];
        st_nZ = c2D.xStart[2];

        loc_gnX = loc_nX + 2;
        loc_gnY = loc_nY + 2;
        loc_gnZ = loc_nZ + 2;

        grid_ndim = loc_nX * loc_nY * loc_nZ;
        grid_gndim = loc_gnX * loc_gnY * loc_gnZ;

        phy_nX = loc_nX - hasPeriodicLayerX + 1;
        phy_nY = loc_nY - hasPeriodicLayerY + 1;
        phy_nZ = loc_nZ - hasPeriodicLayerZ + 1;

        phy_ndim = phy_nX * phy_nY * phy_nZ;

        // GENERATE MPI VECTORS
        MPI_Type_vector(loc_gnY, loc_gnX, loc_gnX, Real_MPI, &XYFace);
        MPI_Type_commit(&XYFace);
        MPI_Type_vector(loc_gnZ, loc_gnX, loc_gnX * loc_gnY, Real_MPI, &XZFace);
        MPI_Type_commit(&XZFace);
        MPI_Type_vector(1, loc_gnX, 0, Real_MPI, &XRow);
        MPI_Type_commit(&XRow);

        MPI_Type_vector(loc_nY, loc_nX, loc_nX, Real_MPI, &XYFace_Periodic);
        MPI_Type_commit(&XYFace_Periodic);
        MPI_Type_vector(loc_nZ, loc_nX, loc_nX * loc_nY, Real_MPI, &XZFace_Periodic);
        MPI_Type_commit(&XZFace_Periodic);

        if (periodicY || !isOnTop) {
            int coord[] = {this_y_pos + 1, this_z_pos};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_north);
        } else { neigh_north = MPI_PROC_NULL; }

        if (periodicY || !isOnBottom) {
            int coord[] = {this_y_pos - 1, this_z_pos};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_south);
        } else { neigh_south = MPI_PROC_NULL; }

        if (periodicZ || !isOnRight) {
            int coord[] = {this_y_pos, this_z_pos + 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_east);
        } else { neigh_east = MPI_PROC_NULL; }

        if (periodicZ || !isOnLeft) {
            int coord[] = {this_y_pos, this_z_pos - 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_west);
        } else { neigh_west = MPI_PROC_NULL; }

        if ((!isOnTop && !isOnRight)
            || (isOnTop && !isOnRight && periodicY)
            || (!isOnTop && isOnRight && periodicZ)) {
            int coord[] = {this_y_pos + 1, this_z_pos + 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_north_east);
        } else { neigh_north_east = MPI_PROC_NULL; }

        if ((!isOnTop && !isOnLeft)
            || (isOnTop && !isOnLeft && periodicY)
            || (!isOnTop && isOnLeft && periodicZ)) {
            int coord[] = {this_y_pos + 1, this_z_pos - 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_north_west);
        } else { neigh_north_west = MPI_PROC_NULL; }

        if ((!isOnBottom && !isOnRight)
            || (isOnBottom && !isOnRight && periodicY)
            || (!isOnBottom && isOnRight && periodicZ)) {
            int coord[] = {this_y_pos - 1, this_z_pos + 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_south_east);
        } else { neigh_south_east = MPI_PROC_NULL; }

        if ((!isOnBottom && !isOnLeft)
            || (isOnBottom && !isOnLeft && periodicY)
            || (!isOnBottom && isOnLeft && periodicZ)) {
            int coord[] = {this_y_pos - 1, this_z_pos - 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_south_west);
        } else { neigh_south_west = MPI_PROC_NULL; }

        // TODO REMOVE, only a check if it is correct
        if (neigh_north != c2D.neighbor[0][2] ||
            neigh_south != c2D.neighbor[0][3] ||
            neigh_east != c2D.neighbor[0][4] ||
            neigh_west != c2D.neighbor[0][5]) {
            std::cerr << "Neighbours are not the same of c2d" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    ~Parameters() {
        MPI_Type_free(&XYFace);
        MPI_Type_free(&XZFace);
        MPI_Type_free(&XRow);
    }
};

#endif //PARAMETERS_HPP
