#ifndef AEROHPC_A_MPIBOUNDARIES_H
#define AEROHPC_A_MPIBOUNDARIES_H

#include "MPITraits.hpp"
#include "Boundaries.hpp"

class MPIBoundaries : public Boundaries {


    std::vector<MPICondition *> _pre_cs;

public:

    MPIBoundaries() = default;

    explicit MPIBoundaries(std::vector<Condition *> cs, std::vector<MPICondition *> pre_cs)
            : Boundaries(std::move(cs)),
              _pre_cs(std::move(pre_cs)) {}

    MPIBoundaries &addPreCond(MPICondition &cond) { _pre_cs.push_back(&cond); return *this; }

    void apply(GridData &grid, Real time) override {
        // Fill all outgoing buffers
        for (MPICondition *pre_c: _pre_cs) pre_c->init(grid);

        // Start exchange of input and output buffers
        for (MPICondition *pre_c: _pre_cs) pre_c->exchange();

        // Apply all boundary condition
        // first we apply all the physical one (this optimizes communication while working)
        for (Condition *bc: _cs) bc->apply(grid, time);

        // then all logical ones
        for (MPICondition *pre_c: _pre_cs) pre_c->apply(grid, time);
    }

};

#endif //AEROHPC_A_MPIBOUNDARIES_H
