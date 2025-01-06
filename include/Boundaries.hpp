#ifndef AEROHPC_A_BOUNDARIES_HPP
#define AEROHPC_A_BOUNDARIES_HPP

#include "Traits.hpp"
#include "Condition.hpp"

class Boundaries {

protected:
    /**
     * Inlet velocity boundary condition
     */
    std::vector<Condition *> _cs;

public:

    Boundaries() = default;

    explicit Boundaries(std::vector<Condition *> cs) : _cs(std::move(cs)) {}

    /**
      * Add a boundary condition to the list of the model
      */
    Boundaries &addCond(Condition &cond) { _cs.push_back(&cond); return *this; }

    /**
      * Apply all the boundary conditions of the model,
      * Be aware that the insertion order of the conditions can affect the result of the application
      */
    virtual void apply(GridData &grid, Real time) { for (Condition *bc: _cs) bc->apply(grid, time); }
};

#endif //AEROHPC_A_BOUNDARIES_HPP
