#ifndef AEROHPC_A_BOUNDARIES_HPP
#define AEROHPC_A_BOUNDARIES_HPP

#include "Traits.hpp"
#include "data/SolverData.hpp"
#include "utils/macroUtils.hpp"

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

#define boundary_add_request_periodic(collection, i, j, k, type_ptr, neigh_rank, face_tag)\
    collection.emplace_back(&rhs_P(i, j, k), type_ptr, neigh_rank, face_tag | P_BUFFER_TAG);



void inline apply_boundaries(Real *data, const Real currentTime, int type) {
    const bool apply_velocity = type & VELOCITY;
    const bool apply_pressure = type & PRESSURE;

    std::vector<MPI_BCRequest> bcs_send;
    std::vector<MPI_BCRequest> bcs_recv;

    // NORTH
    {
        const index_t north_inner_j = params.loc_nY - params.hasPeriodicLayerY - 1;
        const index_t north_ghost_j = params.loc_nY - params.hasPeriodicLayerY;

        const Real Y = real(north_ghost_j + params.st_nY) * params.dY + params.originY;
        // VELOCITY
        if (apply_velocity) {
            if (params.neigh_north != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, north_inner_j, -1, &params.XZFace, params.neigh_north, NORTH_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, north_ghost_j, -1, &params.XZFace, params.neigh_north, SOUTH_BUFFER_TAG);
            } else {
                if (domData.northType == DIRICHLET) {
                    // apply on face
                    ITERATE_XZ_FACE(i, k, true)
                        // On y = phy_dim for domain point we have exact for V
                        V(data, i, north_inner_j, k) = domData.northBF.VF(X + params.dX2, Y, Z + params.dZ2, currentTime);
                        // For ghost points we have useless V, other approximate

                        U(data, i, north_ghost_j, k) = 2 * domData.northBF.UF(X + params.dX, Y, Z + params.dZ2, currentTime)
                                                       - U(data, i, north_inner_j, k);
                        V(data, i, north_ghost_j, k) = 0;
                        W(data, i, north_ghost_j, k) = 2 * domData.northBF.WF(X + params.dX2, Y, Z + params.dZ, currentTime)
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
                boundary_add_request_pressure(bcs_send, data, -1, north_inner_j, -1, &params.XZFace, params.neigh_north, NORTH_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, north_ghost_j, -1, &params.XZFace, params.neigh_north, SOUTH_BUFFER_TAG);
            }
        }
    }

    // SOUTH
    {
        constexpr index_t south_inner_j = 0;
        constexpr index_t south_ghost_j = -1;

        const Real Y = real(south_inner_j + params.st_nY) * params.dY + params.originY;
        // VELOCITY
        if (apply_velocity) {
            if (params.neigh_south != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, south_inner_j, -1, &params.XZFace, params.neigh_south, SOUTH_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, south_ghost_j, -1, &params.XZFace, params.neigh_south, NORTH_BUFFER_TAG);
            } else {
                if (domData.southType == DIRICHLET) {
                    // apply on face
                    ITERATE_XZ_FACE(i, k, true)
                        // On y = 0 for ghost point we hae exact for V, other approximate
                        U(data, i, south_ghost_j, k) = 2 * domData.southBF.UF(X + params.dX, Y, Z + params.dZ2, currentTime)
                                                       - U(data, i, south_inner_j, k);
                        V(data, i, south_ghost_j, k) = domData.southBF.VF(X + params.dX, Y, Z + params.dZ2, currentTime);
                        W(data, i, south_ghost_j, k) = 2 * domData.southBF.WF(X + params.dX2, Y, Z + params.dZ, currentTime)
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
                boundary_add_request_pressure(bcs_send, data, -1, south_inner_j, -1, &params.XZFace, params.neigh_south, SOUTH_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, south_ghost_j, -1, &params.XZFace, params.neigh_south, NORTH_BUFFER_TAG);
            }
        }
    }

    // EAST
    {
        const index_t east_inner_k = params.loc_nZ - params.hasPeriodicLayerZ - 1;
        const index_t east_ghost_k = params.loc_nZ - params.hasPeriodicLayerZ;

        const Real Z = real(east_ghost_k + params.st_nZ) * params.dZ + params.originZ;
        // VELOCITY
        if (apply_velocity) {
            if (params.neigh_east != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, -1, east_inner_k, &params.XYFace, params.neigh_east, EAST_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, -1, east_ghost_k, &params.XYFace, params.neigh_east, WEST_BUFFER_TAG);
            } else {
                if (domData.eastType == DIRICHLET) {
                    // apply on face
                    ITERATE_XY_FACE(i, j, true)
                        // On y = phy_dim for domain point we have exact for V
                        W(data, i, j, east_inner_k) = domData.eastBF.WF(X + params.dX2, Y + params.dY2, Z, currentTime);
                        // For ghost points we have useless V, other approximate

                        U(data, i, j, east_ghost_k) = 2 * domData.eastBF.UF(X + params.dX, Y + params.dY2, Z, currentTime)
                                                      - U(data, i, j, east_inner_k);
                        V(data, i, j, east_ghost_k) = 2 * domData.eastBF.VF(X + params.dX2, Y + params.dY, Z, currentTime)
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
                boundary_add_request_pressure(bcs_send, data, -1, -1, east_inner_k, &params.XYFace, params.neigh_east, EAST_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, -1, east_ghost_k, &params.XYFace, params.neigh_east, WEST_BUFFER_TAG);
            }
        }
    }

    // WEST
    {
        constexpr index_t west_inner_k = 0;
        constexpr index_t west_ghost_k = -1;

        const Real Z = real(west_inner_k + params.st_nZ) * params.dZ + params.originZ;
        // VELOCITY
        if (apply_velocity) {
            if (params.neigh_west != MPI_PROC_NULL) {
                boundary_add_request_velocity(bcs_send, data, -1, -1, west_inner_k, &params.XYFace, params.neigh_west, WEST_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, -1, west_ghost_k, &params.XYFace, params.neigh_west, EAST_BUFFER_TAG);
            } else {
                if (domData.westType == DIRICHLET) {
                    // apply on face
                    ITERATE_XY_FACE(i, j, true)
                        // On y = 0 for ghost point we hae exact for V, other approximate
                        U(data, i, j, west_ghost_k) = 2 * domData.westBF.UF(X + params.dX, Y + params.dY2, Z, currentTime)
                                                      - U(data, i, j, west_inner_k);
                        V(data, i, j, west_ghost_k) = 2 * domData.westBF.VF(X + params.dX2, Y + params.dY, Z, currentTime)
                                                      - V(data, i, j, west_inner_k);
                        W(data, i, j, west_ghost_k) = domData.westBF.WF(X + params.dX2, Y + params.dY2, Z, currentTime);
                    ITERATE_FACE_END()
                } else {
                    // TODO Write neumann condition
                }
            }
        }
        // PRESSURE
        if (apply_pressure) {
            if (params.neigh_west != MPI_PROC_NULL) {
                boundary_add_request_pressure(bcs_send, data, -1, -1, west_inner_k, &params.XYFace, params.neigh_west, WEST_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, -1, west_ghost_k, &params.XYFace, params.neigh_west, EAST_BUFFER_TAG);
            }
        }
    }

    // BACK
    {
        const index_t back_inner_i = params.loc_nX - params.hasPeriodicLayerX - 1;
        const index_t back_ghost_i = params.loc_nX - params.hasPeriodicLayerX;

        constexpr index_t back_periodic_i = 0;

        const Real X = real(back_ghost_i + params.st_nX) * params.dX + params.originX;
        // VELOCITY
        if (apply_velocity) {
            if (params.periodicX) {
                ITERATE_YZ_FACE(j, k, false)
                    U(data, back_ghost_i, j, k) = U(data, back_periodic_i, j, k);
                    V(data, back_ghost_i, j, k) = V(data, back_periodic_i, j, k);
                    W(data, back_ghost_i, j, k) = W(data, back_periodic_i, j, k);
                ITERATE_FACE_END()
            } else {
                if (domData.backType == DIRICHLET) {
                    // apply on face
                    ITERATE_YZ_FACE(j, k, true)
                        // On y = phy_dim for domain point we have exact for V
                        U(data, back_inner_i, j, k) = domData.backBF.UF(X, Y + params.dY2, Z + params.dZ2, currentTime);
                        // For ghost points we have useless V, other approximate
                        U(data, back_ghost_i, j, k) = 0;
                        V(data, back_ghost_i, j, k) = 2 * domData.backBF.VF(X, Y + params.dY, Z + params.dZ2, currentTime)
                                                      - V(data, back_inner_i, j, k);
                        W(data, back_ghost_i, j, k) = 2 * domData.backBF.WF(X, Y + params.dY2, Z + params.dZ, currentTime)
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

        const Real X = real(front_inner_i + params.st_nX) * params.dX + params.originX;
        // VELOCITY
        if (apply_velocity) {
            if (params.periodicX) {
                ITERATE_YZ_FACE(j, k, false)
                    U(data, front_ghost_i, j, k) = U(data, front_periodic_i, j, k);
                    V(data, front_ghost_i, j, k) = V(data, front_periodic_i, j, k);
                    W(data, front_ghost_i, j, k) = W(data, front_periodic_i, j, k);
                ITERATE_FACE_END()
            } else {
                if (domData.frontType == DIRICHLET) {
                    // apply on face
                    ITERATE_YZ_FACE(j, k, true)
                        // On y = 0 for ghost point we hae exact for V, other approximate
                        U(data, front_ghost_i, j, k) = domData.frontBF.UF(X, Y + params.dY2, Z + params.dZ2, currentTime);
                        V(data, front_ghost_i, j, k) = 2 * domData.frontBF.VF(X, Y + params.dY, Z + params.dZ2, currentTime)
                                                       - V(data, front_inner_i, j, k);
                        W(data, front_ghost_i, j, k) = 2 * domData.frontBF.WF(X, Y + params.dY2, Z + params.dZ, currentTime)
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
        const index_t ne_inner_j = params.loc_nY - params.hasPeriodicLayerY - 1;
        const index_t ne_inner_k = params.loc_nZ - params.hasPeriodicLayerZ - 1;
        const index_t ne_ghost_j = params.loc_nY - params.hasPeriodicLayerY;
        const index_t ne_ghost_k = params.loc_nZ - params.hasPeriodicLayerZ;

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
        const index_t ne_inner_j = params.loc_nY - params.hasPeriodicLayerY - 1;
        constexpr index_t ne_inner_k = 0;
        const index_t ne_ghost_j = params.loc_nY - params.hasPeriodicLayerY;
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
        const index_t ne_inner_k = params.loc_nZ - params.hasPeriodicLayerZ - 1;
        constexpr index_t ne_ghost_j = -1;
        const index_t ne_ghost_k = params.loc_nZ - params.hasPeriodicLayerZ;

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
                boundary_add_request_velocity(bcs_send, data, -1, ne_inner_j, ne_inner_k, &params.XRow, params.neigh_south_west, SOUTH_WEST_BUFFER_TAG);
                boundary_add_request_velocity(bcs_recv, data, -1, ne_ghost_j, ne_ghost_k, &params.XRow, params.neigh_south_west, NORTH_EAST_BUFFER_TAG);
            }
        }
        if (apply_pressure) {
            if (params.neigh_south_west != MPI_PROC_NULL) {
                boundary_add_request_pressure(bcs_send, data, -1, ne_inner_j, ne_inner_k, &params.XRow, params.neigh_south_west, SOUTH_WEST_BUFFER_TAG);
                boundary_add_request_pressure(bcs_recv, data, -1, ne_ghost_j, ne_ghost_k, &params.XRow, params.neigh_south_west, NORTH_EAST_BUFFER_TAG);
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

void inline exchange_periodicLayer() {
    std::vector<MPI_BCRequest> bcs_send;
    std::vector<MPI_BCRequest> bcs_recv;

    // NORTH
    {
        const index_t periodic_j = params.loc_nY - 1;

        // PRESSURE
        if (params.periodicY && params.isOnTop) {
            boundary_add_request_periodic(bcs_recv, 0, periodic_j, 0, &params.XZFace_Periodic, params.neigh_north, SOUTH_BUFFER_TAG);
        }
    }

    // SOUTH
    {
        constexpr index_t inner_j = 0;

        // PRESSURE
        if (params.periodicY && params.isOnBottom) {
            boundary_add_request_periodic(bcs_send, 0, inner_j, 0, &params.XZFace_Periodic, params.neigh_south, SOUTH_BUFFER_TAG);
        }
    }

    // EAST
    {
        const index_t periodic_k = params.loc_nZ - 1;

        // PRESSURE
        if (params.periodicZ && params.isOnRight) {
            boundary_add_request_periodic(bcs_recv, 0, 0, periodic_k, &params.XYFace_Periodic, params.neigh_east, WEST_BUFFER_TAG);
        }
    }

    // WEST
    {
        constexpr index_t inner_k = 0;

        // PRESSURE
        if (params.periodicZ && params.isOnLeft) {
            boundary_add_request_periodic(bcs_send, 0, 0, inner_k, &params.XYFace_Periodic, params.neigh_west, WEST_BUFFER_TAG);
        }
    }

    // FRONT
    {
        constexpr index_t inner_i = 0;

        const index_t periodic_i = params.loc_nX - 1;

        if (params.periodicX) {
            ITERATE_YZ_FACE(j, k, false)
                rhs_P(periodic_i, j, k) = rhs_P(inner_i, j, k);
            ITERATE_FACE_END()
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
}

#endif //AEROHPC_A_BOUNDARIES_HPP
