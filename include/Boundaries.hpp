#ifndef AEROHPC_A_BOUNDARIES_HPP
#define AEROHPC_A_BOUNDARIES_HPP

#include "Traits.hpp"
#include "data/SolverData.hpp"
#include "utils/macroUtils.hpp"

#define TYPE_VELOCITY 1
#define TYPE_PRESSURE 2

class MPI_BCRequest {
public:
    Real *data_basePtr;
    const MPI_Datatype *datatype;
    const int neighRank;
    const int bufferTag;

    MPI_Request request{};

    MPI_BCRequest(Real *data_basePtr, const MPI_Datatype *datatype, const int neighRank, const int bufferTag)
        : data_basePtr(data_basePtr), datatype(datatype), neighRank(neighRank), bufferTag(bufferTag) {
    }
};

#define boundary_add_request_velocity(collection, data_ptr, i, j, k, type_ptr, neigh_rank, face_tag)\
    collection.emplace_back(&U(data_ptr, i, j, k), type_ptr, neigh_rank, face_tag | U_BUFFER_TAG);\
    collection.emplace_back(&V(data_ptr, i, j, k), type_ptr, neigh_rank, face_tag | V_BUFFER_TAG);\
    collection.emplace_back(&W(data_ptr, i, j, k), type_ptr, neigh_rank, face_tag | W_BUFFER_TAG)

#define boundary_add_request_pressure(collection, data_ptr, i, j, k, type_ptr, neigh_rank, face_tag)\
    collection.emplace_back(&P(data_ptr, i, j, k), type_ptr, neigh_rank, face_tag | P_BUFFER_TAG);

void inline apply_boundaries(Real *data, const Real currentTime, int type) {
    const bool apply_velocity = type & TYPE_VELOCITY;
    const bool apply_pressure = type & TYPE_PRESSURE;

    std::vector<MPI_BCRequest> bcs_send;
    std::vector<MPI_BCRequest> bcs_recv;

    // NORTH
    {
        const index_t north_inner_j = params.loc_nY - 1;
        const index_t north_ghost_j = params.loc_nY;

        const Real y = real(north_ghost_j + params.st_nY) * params.dY + params.originY;
        // VELOCITY
        if (apply_velocity) {
            if (params.neigh_north != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, north_inner_j, 0, &params.XZFace, params.neigh_north, NORTH_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, north_ghost_j, 0, &params.XZFace, params.neigh_north, SOUTH_BUFFER_TAG);
            } else {
                if (domData.northType == BOUNDARY_DIRICHLET) {
                    // apply on face
                    ITERATE_XZ_FACE(i, k, x, z)
                        // On y = phy_dim for domain point we have exact for V
                        V(data, i, north_inner_j, k) = domData.northBF.VF(x + params.dX2, y, z + params.dZ2, currentTime);
                        // For ghost points we have useless V, other approximate

                        U(data, i, north_ghost_j, k) = 2 * domData.northBF.UF(x + params.dX, y, z + params.dZ2, currentTime)
                                                       - U(data, i, north_inner_j, k);
                        V(data, i, north_ghost_j, k) = 0;
                        W(data, i, north_ghost_j, k) = 2 * domData.northBF.WF(x + params.dX2, y, z + params.dZ, currentTime)
                                                       - W(data, i, north_inner_j, k);
                    ITERATE_FACE_END()
                } else {
                    // TODO Write neumann condition
                }
            }
        }
        // PRESSURE
        if (apply_pressure) {
            if (params.neigh_north != MPI_PROC_NULL) {
                boundary_add_request_pressure(bcs_send, data, -1, north_inner_j, 0, &params.XZFace, params.neigh_north, NORTH_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, north_ghost_j, 0, &params.XZFace, params.neigh_north, SOUTH_BUFFER_TAG);
            }
        }
    }

    // SOUTH
    {
        constexpr index_t south_inner_j = 0;
        constexpr index_t south_ghost_j = -1;

        const Real y = real(south_inner_j + params.st_nY) * params.dY + params.originY;
        // VELOCITY
        if (apply_velocity) {
            if (params.neigh_south != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, south_inner_j, 0, &params.XZFace, params.neigh_south, SOUTH_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, south_ghost_j, 0, &params.XZFace, params.neigh_south, NORTH_BUFFER_TAG);
            } else {
                if (domData.southType == BOUNDARY_DIRICHLET) {
                    // apply on face
                    ITERATE_XZ_FACE(i, k, x, z)
                        // On y = 0 for ghost point we hae exact for V, other approximate
                        U(data, i, south_ghost_j, k) = 2 * domData.southBF.UF(x + params.dX, y, z + params.dZ2, currentTime)
                                                       - U(data, i, south_inner_j, k);
                        V(data, i, south_ghost_j, k) = domData.southBF.VF(x + params.dX, y, z + params.dZ2, currentTime);
                        W(data, i, south_ghost_j, k) = 2 * domData.southBF.WF(x + params.dX2, y, z + params.dZ, currentTime)
                                                       - W(data, i, south_inner_j, k);
                    ITERATE_FACE_END()
                } else {
                    // TODO Write neumann condition
                }
            }
        }
        // PRESSURE
        if (apply_pressure) {
            if (params.neigh_south != MPI_PROC_NULL) {
                boundary_add_request_pressure(bcs_send, data, -1, south_inner_j, 0, &params.XZFace, params.neigh_south, SOUTH_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, south_ghost_j, 0, &params.XZFace, params.neigh_south, NORTH_BUFFER_TAG);
            }
        }
    }

    // EAST
    {
        const index_t east_inner_k = params.loc_nZ - 1;
        const index_t east_ghost_k = params.loc_nZ;

        const Real z = real(east_ghost_k + params.st_nZ) * params.dZ + params.originZ;
        // VELOCITY
        if (apply_velocity) {
            if (params.neigh_east != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, 0, east_inner_k, &params.XYFace, params.neigh_east, EAST_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, 0, east_ghost_k, &params.XYFace, params.neigh_east, WEST_BUFFER_TAG);
            } else {
                if (domData.eastType == BOUNDARY_DIRICHLET) {
                    // apply on face
                    ITERATE_XY_FACE(i, j, x, y)
                        // On y = phy_dim for domain point we have exact for V
                        W(data, i, j, east_inner_k) = domData.eastBF.WF(x + params.dX2, y + params.dY2, z, currentTime);
                        // For ghost points we have useless V, other approximate

                        U(data, i, j, east_ghost_k) = 2 * domData.eastBF.UF(x + params.dX, y + params.dY2, z, currentTime)
                                                      - U(data, i, j, east_inner_k);
                        V(data, i, j, east_ghost_k) = 2 * domData.eastBF.VF(x + params.dX2, y + params.dY, z, currentTime)
                                                      - V(data, i, j, east_inner_k);
                        W(data, i, j, east_ghost_k) = 0;
                    ITERATE_FACE_END()
                } else {
                    // TODO Write neumann condition
                }
            }
        }
        // PRESSURE
        if (apply_pressure) {
            if (params.neigh_east != MPI_PROC_NULL) {
                boundary_add_request_pressure(bcs_send, data, -1, 0, east_inner_k, &params.XYFace, params.neigh_east, EAST_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, 0, east_ghost_k, &params.XYFace, params.neigh_east, WEST_BUFFER_TAG);
            }
        }
    }

    // WEST
    {
        constexpr index_t west_inner_k = 0;
        constexpr index_t west_ghost_k = -1;

        const Real z = real(west_inner_k + params.st_nZ) * params.dZ + params.originZ;
        // VELOCITY
        if (apply_velocity) {
            if (params.neigh_west != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, 0, west_inner_k, &params.XYFace, params.neigh_west, WEST_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, 0, west_ghost_k, &params.XYFace, params.neigh_west, EAST_BUFFER_TAG);
            } else {
                if (domData.westType == BOUNDARY_DIRICHLET) {
                    // apply on face
                    ITERATE_XY_FACE(i, j, x, y)
                        // On y = 0 for ghost point we hae exact for V, other approximate
                        U(data, i, j, west_ghost_k) = 2 * domData.westBF.UF(x + params.dX, y + params.dY2, z, currentTime)
                                                      - U(data, i, j, west_inner_k);
                        V(data, i, j, west_ghost_k) = 2 * domData.westBF.VF(x + params.dX2, y + params.dY, z, currentTime)
                                                      - V(data, i, j, west_inner_k);
                        W(data, i, j, west_ghost_k) = domData.westBF.WF(x + params.dX2, y + params.dY2, z, currentTime);
                    ITERATE_FACE_END()
                } else {
                    // TODO Write neumann condition
                }
            }
        }
        // PRESSURE
        if (apply_pressure) {
            if (params.neigh_west != MPI_PROC_NULL) {
                boundary_add_request_pressure(bcs_send, data, -1, 0, west_inner_k, &params.XYFace, params.neigh_west, WEST_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, 0, west_ghost_k, &params.XYFace, params.neigh_west, EAST_BUFFER_TAG);
            }
        }
    }

    // BACK
    {
        const index_t back_inner_i = params.loc_nX - 1;
        const index_t back_ghost_i = params.loc_nX;

        constexpr index_t back_periodic_i = 0;

        const Real x = real(back_ghost_i + params.st_nX) * params.dX + params.originX;
        // VELOCITY
        if (apply_velocity) {
            if (params.periodicX) {
                ITERATE_YZ_FACE(j, k, y, z)
                    U(data, back_ghost_i, j, k) = U(data, back_periodic_i, j, k);
                    V(data, back_ghost_i, j, k) = V(data, back_periodic_i, j, k);
                    W(data, back_ghost_i, j, k) = W(data, back_periodic_i, j, k);
                ITERATE_FACE_END()
            } else {
                if (domData.backType == BOUNDARY_DIRICHLET) {
                    // apply on face
                    ITERATE_YZ_FACE(j, k, y, z)
                        // On y = phy_dim for domain point we have exact for V
                        U(data, back_inner_i, j, k) = domData.backBF.UF(x, y + params.dY2, z + params.dZ2, currentTime);
                        // For ghost points we have useless V, other approximate
                        U(data, back_ghost_i, j, k) = 0;
                        V(data, back_ghost_i, j, k) = 2 * domData.backBF.VF(x, y + params.dY, z + params.dZ2, currentTime)
                                                      - V(data, back_inner_i, j, k);
                        W(data, back_ghost_i, j, k) = 2 * domData.backBF.WF(x, y + params.dY2, z + params.dZ, currentTime)
                                                      - W(data, back_inner_i, j, k);
                    ITERATE_FACE_END()
                } else {
                    // TODO Write neumann condition
                }
            }
        }
    }

    // FRONT
    {
        constexpr index_t front_inner_i = 0;
        constexpr index_t front_ghost_i = -1;

        const index_t front_periodic_i = params.loc_nX - 1;

        const Real x = real(front_inner_i + params.st_nX) * params.dX + params.originX;
        // VELOCITY
        if (apply_velocity) {
            if (params.periodicX) {
                ITERATE_YZ_FACE(j, k, y, z)
                    U(data, front_ghost_i, j, k) = U(data, front_periodic_i, j, k);
                    V(data, front_ghost_i, j, k) = V(data, front_periodic_i, j, k);
                    W(data, front_ghost_i, j, k) = W(data, front_periodic_i, j, k);
                ITERATE_FACE_END()
            } else {
                if (domData.frontType == BOUNDARY_DIRICHLET) {
                    // apply on face
                    ITERATE_YZ_FACE(j, k, y, z)
                        // On y = 0 for ghost point we hae exact for V, other approximate
                        U(data, front_ghost_i, j, k) = domData.frontBF.UF(x, y + params.dY2, z + params.dZ2, currentTime);
                        V(data, front_ghost_i, j, k) = 2 * domData.frontBF.VF(x, y + params.dY, z + params.dZ2, currentTime)
                                                       - V(data, front_inner_i, j, k);
                        W(data, front_ghost_i, j, k) = 2 * domData.frontBF.WF(x, y + params.dY2, z + params.dZ, currentTime)
                                                       - W(data, front_inner_i, j, k);
                    ITERATE_FACE_END()
                } else {
                    // TODO Write neumann condition
                }
            }
        }
    }


    // NORTH EAST
    {
        const index_t ne_inner_j = params.loc_nY - 1;
        const index_t ne_inner_k = params.loc_nZ - 1;
        const index_t ne_ghost_j = params.loc_nY;
        const index_t ne_ghost_k = params.loc_nZ;

        if (apply_velocity) {
            if (params.neigh_north_east != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, ne_inner_j, ne_inner_k, &params.XRow, params.neigh_north_east, NORTH_EAST_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, ne_ghost_j, ne_ghost_k, &params.XRow, params.neigh_north_east, SOUTH_WEST_BUFFER_TAG);
            }
        }
        if (apply_pressure) {
            if (params.neigh_north_east != MPI_PROC_NULL) {
                boundary_add_request_pressure(bcs_send, data, -1, ne_inner_j, ne_inner_k, &params.XRow, params.neigh_north_east, NORTH_EAST_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, ne_ghost_j, ne_ghost_k, &params.XRow, params.neigh_north_east, SOUTH_WEST_BUFFER_TAG);
            }
        }
    }

    // NORTH WEST
    {
        const index_t ne_inner_j = params.loc_nY - 1;
        constexpr index_t ne_inner_k = 0;
        const index_t ne_ghost_j = params.loc_nY;
        constexpr index_t ne_ghost_k = -1;

        if (apply_velocity) {
            if (params.neigh_north_west != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, ne_inner_j, ne_inner_k, &params.XRow, params.neigh_north_west, NORTH_WEST_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, ne_ghost_j, ne_ghost_k, &params.XRow, params.neigh_north_west, SOUTH_EAST_BUFFER_TAG);
            }
        }
        if (apply_pressure) {
            if (params.neigh_north_west != MPI_PROC_NULL) {
                boundary_add_request_pressure(bcs_send, data, -1, ne_inner_j, ne_inner_k, &params.XRow, params.neigh_north_west, NORTH_WEST_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, ne_ghost_j, ne_ghost_k, &params.XRow, params.neigh_north_west, SOUTH_EAST_BUFFER_TAG);
            }
        }
    }


    // SOUTH EAST
    {
        constexpr index_t ne_inner_j = 0;
        const index_t ne_inner_k = params.loc_nZ - 1;
        constexpr index_t ne_ghost_j = -1;
        const index_t ne_ghost_k = params.loc_nZ;

        if (apply_velocity) {
            if (params.neigh_south_east != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, ne_inner_j, ne_inner_k, &params.XRow, params.neigh_south_east, SOUTH_EAST_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, ne_ghost_j, ne_ghost_k, &params.XRow, params.neigh_south_east, NORTH_WEST_BUFFER_TAG);
            }
        }
        if (apply_pressure) {
            if (params.neigh_south_east != MPI_PROC_NULL) {
                boundary_add_request_pressure(bcs_send, data, -1, ne_inner_j, ne_inner_k, &params.XRow, params.neigh_south_east, SOUTH_EAST_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, ne_ghost_j, ne_ghost_k, &params.XRow, params.neigh_south_east, NORTH_WEST_BUFFER_TAG);
            }
        }
    }

    // SOUTH WEST
    {
        constexpr index_t ne_inner_j = 0;
        constexpr index_t ne_inner_k = 0;
        constexpr index_t ne_ghost_j = -1;
        constexpr index_t ne_ghost_k = -1;

        if (apply_velocity) {
            if (params.neigh_south_west != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, ne_inner_j, ne_inner_k, &params.XRow, params.neigh_south_east, SOUTH_WEST_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, ne_ghost_j, ne_ghost_k, &params.XRow, params.neigh_south_east, NORTH_EAST_BUFFER_TAG);
            }
        }
        if (apply_pressure) {
            if (params.neigh_south_west != MPI_PROC_NULL) {
                boundary_add_request_pressure(bcs_send, data, -1, ne_inner_j, ne_inner_k, &params.XRow, params.neigh_south_east, SOUTH_WEST_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, ne_ghost_j, ne_ghost_k, &params.XRow, params.neigh_south_east, NORTH_EAST_BUFFER_TAG);
            }
        }
    }

    // EXCHANGE FACES WITH OTHER SUBDOMAINS
    for (MPI_BCRequest &req: bcs_send) {
        MPI_Isend(req.data_basePtr, 1, *req.datatype, req.neighRank, req.bufferTag, MPI_COMM_WORLD, &req.request);
    }
    for (MPI_BCRequest &req: bcs_recv) {
        MPI_Irecv(req.data_basePtr, 1, *req.datatype, req.neighRank, req.bufferTag, MPI_COMM_WORLD, &req.request);
    }
    for (MPI_BCRequest &req: bcs_send) {
        MPI_Wait(&req.request, MPI_STATUS_IGNORE);
    }
    for (MPI_BCRequest &req: bcs_recv) {
        MPI_Wait(&req.request, MPI_STATUS_IGNORE);
    }
}


#endif //AEROHPC_A_BOUNDARIES_HPP
