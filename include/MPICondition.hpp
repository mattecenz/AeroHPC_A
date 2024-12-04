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
    typedef void (*BufferExchanger)(GridData &bufferOut, GridData &bufferIn, MPI_Request *requestOut, MPI_Request *requestIn, int proc_rank);

    /**
     * Typedef shortening lambda definition of mapping function
     */
    typedef void (*BufferMapper)(GridData &grid, GridData &buffer, MPI_Request *requestIn, MPI_Request *requestOut);


    MPICondition() = delete;


    MPICondition(const BufferInitializer &initializer,
                 const BufferExchanger &exchanger,
                 const BufferMapper &mapper,
                 const GridStructure &bufferStructure,
                 const int proc_rank) : _initializer(initializer),
                                        _exchanger(exchanger),
                                        _mapper(mapper),
                                        _proc_rank(proc_rank),
                                        _bufferIn(bufferStructure, false),
                                        _bufferOut(bufferStructure, false),
                                        _requestIn(),
                                        _requestOut() {}


    /**
     *
     */
    void init(GridData &grid) { _initializer(grid, _bufferOut); }

    void exchange() { _exchanger(_bufferOut, _bufferIn, &_requestOut, &_requestIn, _proc_rank); }

    void apply(GridData &grid, const Real time) override { _mapper(grid, _bufferOut, &_requestIn, &_requestOut); }

private:

    const BufferMapper _mapper;
    const BufferInitializer _initializer;
    const BufferExchanger _exchanger;

    const int _proc_rank;

    MPI_Request _requestIn, _requestOut;

    GridData _bufferIn, _bufferOut;
};

#endif //AEROHPC_A_MPICONDITION_H
