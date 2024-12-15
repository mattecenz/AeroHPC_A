#ifndef AEROHPC_A_MPICONDITION_H
#define AEROHPC_A_MPICONDITION_H

#include "MPITraits.hpp"
#include "Condition.hpp"

class MPICondition : public Condition {
public:

    /**
     * Typedef shortening lambda definition of MPI buffer initializer
     */
    typedef void (*BufferInitializer)(GridData &grid, GridData &bufferOut);

    /**
     * Typedef shortening lambda definition of MPI buffer exchanger
     */
    typedef void (*BufferExchanger)(GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut, MPI_Request *requestIn, int neigh_rank);

    /**
     * Typedef shortening lambda definition of mapping function
     */
    typedef void (*BufferMapper)(GridData &grid, GridData &buffer, MPI_Request *requestOut, MPI_Request *requestIn);


    MPICondition() = delete;

    /**
     * Construct an MPI Boundary condition with given Buffer manages,
     * setup ingoing and outgoing buffers with given structure,
     * store information of neighbour rank
     */
    MPICondition(const BufferInitializer &initializer,
                 const BufferExchanger &exchanger,
                 const BufferMapper &mapper,
                 const GridStructure &bufferStructure,
                 const int neigh_rank) : _initializer(initializer),
                                         _exchanger(exchanger),
                                         _mapper(mapper),
                                         _neigh_rank(neigh_rank),
                                         _bufferIn(bufferStructure, false),
                                         _bufferOut(bufferStructure, false){}


    /**
     * This function call the buffer initializer, this will setup outgoing buffer with given grid data
     */
    void init(GridData &grid) { _initializer(grid, _bufferOut); }

    /**
     * This function call the buffer exchanger, start MPI communications
     */
    void exchange() { _exchanger(_bufferOut, _bufferIn, &_requestOut, &_requestIn, _neigh_rank); }

    /**
     * This function call the grid mapper, this will load ingoing buffer data into the given grid
     */
    void apply(GridData &grid, const Real time) override { _mapper(grid, _bufferIn, &_requestOut, &_requestIn); }

private:

    /**
     * The grid mapper lambda
     */
    const BufferMapper _mapper;
    /**
     * The buffer initializer lambda
     */
    const BufferInitializer _initializer;
    /**
     * The buffer exchanger lambda
     */
    const BufferExchanger _exchanger;

    /**
     * The rank of the neighbour the MPI communications are done with
     */
    const int _neigh_rank;

    /**
     * MPI request objects, they store info about Communication status,
     * they are mainly used as waiting lock for communication
     */
    MPI_Request _requestIn{}, _requestOut{};

    /**
     * The ingoing and outgoing buffer
     */
    GridData _bufferIn, _bufferOut;
};

#endif //AEROHPC_A_MPICONDITION_H
