#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "Traits.hpp"
#include "C2Decomp.hpp"

class Parameters {
public:
    // Spatial Info
    Real dimX, dimY, dimZ;
    Real originX, originY, originZ;

    // Logical Grid Info
    index_t glob_nX, glob_nY, glob_nZ; // Number of gloabl nodes
    index_t loc_nX; index_t loc_nY; index_t loc_nZ; // Number of local nodes
    index_t st_nX; index_t st_nY; index_t st_nZ; // Local nodes starting global index
    index_t loc_gnX; index_t loc_gnY; index_t loc_gnZ; // Number of loocal nodes + ghosts
    index_t grid_ndim; // Total number of local nodes
    index_t grid_gndim; // Total number of local nodes + ghosts

    // Physical Grid Info
    Real dX; Real dY; Real dZ; // Physical spacing between nodes
    Real dX2; Real dY2; Real dZ2; // Physical spacing between staggered nodes

    // Domain Info
    int neigh_front, neigh_back,
            neigh_north, neigh_south,
            neigh_east, neigh_west,
            neigh_north_east, neigh_north_west,
            neigh_south_east, neigh_south_west;

    // Boundary Info
    bool isOnTop, isOnBottom, isOnRight, isOnLeft;
    bool periodicX, periodicY, periodicZ;

    // MPI Types
    MPI_Datatype XYFace;
    MPI_Datatype XZFace;
    MPI_Datatype XRow;

    // Time Info
    Real dt; // Physical delta Time between timesteps

    // Solver Info
    index_t timesteps;
    Real Re;

    Parameters(const Real dimX, const Real dimY, const Real dimZ,
               const Real originX, const Real originY, const Real originZ,
               const Real glob_nX, const Real glob_nY, const Real glob_nZ,
               const Real dt, const index_t timesteps, const Real Re,
                const bool periodicBC[3], const C2Decomp &c2D)
        : dimX(dimX), dimY(dimY), dimZ(dimZ),
          originX(originX), originY(originY), originZ(originZ),
          glob_nX(glob_nX), glob_nY(glob_nY), glob_nZ(glob_nZ),
          dt(dt), timesteps(timesteps), Re(Re) {
        dX = dimX / real(glob_nX);
        dY = dimY / real(glob_nY);
        dZ = dimZ / real(glob_nZ);

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

        // GENERATE MPI VECTORS
        MPI_Type_vector(loc_nY, loc_gnX, loc_gnX, Real_MPI, &XYFace);
        MPI_Type_commit(&XYFace);
        MPI_Type_vector(loc_nZ, loc_gnX, loc_gnX * loc_gnY, Real_MPI, &XYFace);
        MPI_Type_commit(&XZFace);
        MPI_Type_vector(1, loc_gnX, 0, Real_MPI, &XRow);
        MPI_Type_commit(&XRow);

        periodicX = periodicBC[0];
        periodicY = periodicBC[1];
        periodicZ = periodicBC[2];

        // FIND NEIGHBORS
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

        {
            int coord[] = {this_y_pos + 1, this_z_pos};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_north);
        }
        {
            int coord[] = {this_y_pos - 1, this_z_pos};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_south);
        }
        {
            int coord[] = {this_y_pos, this_z_pos + 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_east);
        }
        {
            int coord[] = {this_y_pos, this_z_pos - 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_west);
        }
        {
            int coord[] = {this_y_pos + 1, this_z_pos + 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_north_east);
        }
        {
            int coord[] = {this_y_pos + 1, this_z_pos - 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_north_west);
        }
        {
            int coord[] = {this_y_pos - 1, this_z_pos + 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_south_east);
        }
        {
            int coord[] = {this_y_pos - 1, this_z_pos - 1};
            MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_south_west);
        }

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
