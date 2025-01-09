#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "Traits.hpp"
#include "C2Decomp.hpp"

#define define3d(type, prefix, suffix) \
    type prefix##X##suffix; type prefix##Y##suffix; type prefix##Z##suffix
#define use3d(type, prefix, suffix) \
    type prefix##X##suffix, type prefix##Y##suffix, type prefix##Z##suffix

class Parameters {
public:
    // Spatial Info
    define3d(Real, dim); // Physical dimension of the global space
    define3d(Real, origin); // Physical origin of the global space

    // Logical Grid Info
    define3d(index_t, glob_n); // Number of gloabl nodes
    define3d(index_t, loc_n); // Number of local nodes
    define3d(index_t, st_n); // Local nodes starting global index
    define3d(index_t, loc_gn); // Number of loocal nodes + ghosts
    index_t grid_ndim; // Total number of local nodes
    index_t grid_gndim; // Total number of local nodes + ghosts

    // Physical Grid Info
    define3d(Real, d); // Physical spacing between nodes
    define3d(Real, d, 2); // Physical spacing between staggered nodes

    // Domain Info
    int neigh_front, neigh_back,
            neigh_north, neigh_south,
            neigh_east, neigh_west,
            neigh_north_east, neigh_north_west,
            neigh_south_east, neigh_south_west;

    // Boundary Info


    // Time Info
    Real dt; // Physical delta Time between timesteps

    // Solver Info
    index_t timesteps;
    Real Re;

    Parameters(use3d(const Real, dim),
               use3d(const Real, origin),
               use3d(const Real, glob_n),
               const Real dt, const index_t timesteps, const Real Re,
               const C2Decomp &c2D)
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

        findNeighbours(c2D);
    }

    void findNeighbours(const C2Decomp &c2D) {
        // global position of this process
        const int n_y_proc = c2D.dims[0];
        const int n_z_proc = c2D.dims[1];
        const int this_y_pos = c2D.coord[0];
        const int this_z_pos = c2D.coord[1];

        // flags to define if the processor is on a physical boundary
        const bool isOnTop = (this_y_pos == n_y_proc - 1);
        const bool isOnBottom = (this_y_pos == 0);
        const bool isOnLeft = (this_z_pos == 0);
        const bool isOnRight = (this_z_pos == n_z_proc - 1);

        int coord[] = {this_y_pos + 1, this_z_pos};
        MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_north);

        coord(this_y_pos - 1, this_z_pos);
        MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_south);

        coord(this_y_pos, this_z_pos + 1);
        MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_east);

        coord(this_y_pos, this_z_pos - 1);
        MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_west);

        coord(this_y_pos + 1, this_z_pos + 1);
        MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_north_east);

        coord(this_y_pos + 1, this_z_pos - 1);
        MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_north_west);

        coord(this_y_pos - 1, this_z_pos + 1);
        MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_south_east);

        coord(this_y_pos - 1, this_z_pos - 1);
        MPI_Cart_rank(c2D.DECOMP_2D_COMM_CART_X, coord, &neigh_south_west);

        // TODO REMOVE, only a check if it is correct
        if (neigh_north != c2D.neighbor[0][2] ||
            neigh_south != c2D.neighbor[0][3] ||
            neigh_east != c2D.neighbor[0][4] ||
            neigh_west != c2D.neighbor[0][5]) {
            std::cerr << "Neighbours are not the same of c2d" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }
};

#endif //PARAMETERS_HPP
