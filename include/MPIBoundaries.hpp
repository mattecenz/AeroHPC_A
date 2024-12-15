#ifndef AEROHPC_A_MPIBOUNDARIES_H
#define AEROHPC_A_MPIBOUNDARIES_H

#include "MPITraits.hpp"
#include "Boundaries.hpp"

class MPIBoundaries : public Boundaries {

    /**
     * Collection of boundary function that needs MPI communications
     */
    std::vector<MPICondition *> _mpi_cs;

public:

    MPIBoundaries() = default;

    /**
     * Construct an MPIBoundaries object initialized with given BC and MPI_BC
     */
    MPIBoundaries(std::vector<Condition *> cs, std::vector<MPICondition *> mpi_cs)
            : Boundaries(std::move(cs)),
              _mpi_cs(std::move(mpi_cs)) {}

    /**
     * Add an MPI BC to the collection
     */
    MPIBoundaries &addMPICond(MPICondition &cond) { _mpi_cs.push_back(&cond); return *this; }

    void apply(GridData &grid, Real time) override {

        // Fill all outgoing buffers
        for (MPICondition *mpi_c: _mpi_cs) mpi_c->init(grid);

        // Start exchange of input and output buffers
        for (MPICondition *mpi_c: _mpi_cs) mpi_c->exchange();

        // then all logical ones
        for (MPICondition *mpi_c: _mpi_cs) mpi_c->apply(grid, time);

        // Apply all boundary condition
        // first we apply all the physical one (this optimizes communication while working)
        for (Condition *bc: _cs) bc->apply(grid, time);


        // Await all processor has applied BC
        MPI_Barrier(MPI_COMM_WORLD);
    }

};

#endif //AEROHPC_A_MPIBOUNDARIES_H
