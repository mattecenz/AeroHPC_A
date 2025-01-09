#ifndef AEROHPC_A_BOUNDARIES_HPP
#define AEROHPC_A_BOUNDARIES_HPP

#include "Traits.hpp"
#include "data/SolverData.hpp"

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

#define add_request_velocity(collection, data_ptr, i, j, k, type_ptr, neigh_rank, face_tag)\
    collection.emplace_back(&U(data_ptr, i, j, k), type_ptr, neigh_rank, face_tag | U_BUFFER_TAG);\
    collection.emplace_back(&V(data_ptr, i, j, k), type_ptr, neigh_rank, face_tag | V_BUFFER_TAG);\
    collection.emplace_back(&W(data_ptr, i, j, k), type_ptr, neigh_rank, face_tag | W_BUFFER_TAG)

#define add_request_pressure(collection, data_ptr, i, j, k, type_ptr, neigh_rank, face_tag)\
    collection.emplace_back(&P(data_ptr, i, j, k), type_ptr, neigh_rank, face_tag | P_BUFFER_TAG);

void inline apply_boundaries(Real *data, const Real currentTime, int type) {
    const bool apply_velocity = (type & TYPE_VELOCITY) != 0;
    const bool apply_pressure = (type & TYPE_PRESSURE) != 0;

    std::vector<MPI_BCRequest> bcs_send;
    std::vector<MPI_BCRequest> bcs_recv;

    const index_t north_inner_j = params.loc_nY - 1;
    const index_t north_ghost_j = params.loc_nY;

    // VELOCITY NORTH
    if constexpr (apply_velocity) {
        if (params.periodicY || !params.isOnTop) {
            add_request_velocity(bcs_send, data, 0, north_inner_j, 0, &params.XZFace, params.neigh_north, NORTH_BUFFER_TAG);
            add_request_velocity(bcs_recv, data, 0, north_ghost_j, 0, &params.XZFace, params.neigh_north, NORTH_BUFFER_TAG);
        } else {
            if (bcsFun.northType == BOUNDARY_DIRICHLET) {
                const Real y = real(north_ghost_j + params.st_nX) * params.dY + params.originY;

                // apply on face
                for (index_t k = 0; k < params.loc_nZ; k++) {
                    for (index_t i = 0; i < params.loc_nX; i++) {
                        Real x = real(i + params.st_nX) * params.dX + params.originX;
                        Real z = real(k + params.st_nZ) * params.dZ + params.originZ;

                        // On y = phy_dim for domain point we have exact for V
                        V(data, i, north_inner_j, k) = bcsFun.northFunction(x + params.dX2, y, z + params.dZ2, currentTime);
                        // For ghost points we have useless V, other approximate

                        U(data, i, north_ghost_j, k) = 2 * bcsFun.northFunction(x + params.dX, y, z + params.dZ2, currentTime)
                                                       - U(data, i, north_inner_j, k);
                        V(data, i, north_ghost_j, k) = 0;
                        W(data, i, north_ghost_j, k) = 2 * bcsFun.northFunction(x + params.dX2, y, z + params.dZ, currentTime)
                                                       - W(datai, north_inner_j, k);
                    }
                }
            } else {
                // TODO Write neumann condition
            }
        }
    }
    // PRESSURE NORTH
    if constexpr (apply_pressure) {
        if (params.periodicY || !params.isOnTop) {
            add_request_pressure(bcs_send, data, 0, north_inner_j, 0, &params.XZFace, params.neigh_north, NORTH_BUFFER_TAG);
            add_request_pressure(bcs_recv, data, 0, north_ghost_j, 0, &params.XZFace, params.neigh_north, NORTH_BUFFER_TAG);
        } else {
            // TODO Write neumann condition
        }
    }

    // TODO WRITE ALL OTHER FACES


    // EXCHANGE FACES WITH OTHER SUBDOMAINS
    for (auto req : bcs_send) {
        MPI_Isend(req.data_basePtr, 1, *req.datatype, req.neighRank, req.bufferTag, MPI_COMM_WORLD, &req.request);
    }
    for (auto req : bcs_recv) {
        MPI_Irecv(req.data_basePtr, 1, *req.datatype, req.neighRank, req.bufferTag, MPI_COMM_WORLD, &req.request);
    }
    for (auto req : bcs_send) {
        MPI_Wait(&req.request, MPI_STATUS_IGNORE);
    }
    for (auto req : bcs_recv) {
        MPI_Wait(&req.request, MPI_STATUS_IGNORE);
    }
}


#endif //AEROHPC_A_BOUNDARIES_HPP
