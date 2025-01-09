#ifndef BOUNDARIESDATA_HPP
#define BOUNDARIESDATA_HPP

#include "Traits.hpp"

class MPI_BCRequest {
public:
    Real *data_basePtr;
    MPI_Datatype *datatype;
    MPI_Request *request;
};

class BoundariesData {
public:
    std::vector<MPI_BCRequest> bcs_send;
    std::vector<MPI_BCRequest> bcs_recv;

    BoundariesData() {

    }
};

#endif //BOUNDARIESDATA_HPP
