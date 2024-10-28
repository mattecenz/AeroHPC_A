#ifndef AEROHPC_A_BOUNDARIES_HPP
#define AEROHPC_A_BOUNDARIES_HPP

#include "Traits.hpp"
#include "Condition.hpp"

template<Addressing_T A>
class Boundaries {

private:
    /**
     * Inlet velocity boundary condition
     */
    std::vector<Condition<A> *> _cs;

public:

    Boundaries() = default;

    Boundaries(std::vector<Condition<A> *> cs) : _cs(cs) {}

    /**
      * Add a boundary condition to the list of the model
      */
    void addCond(Condition<A> &cond) { _cs.push_back(&cond); }

    /**
      * Apply all the boundary conditions of the model,
      * Be aware that the insertion order of the conditions can affect the result of the application
      */
    void apply(Grid<A> &grid, Real time) { for (Condition<A> *bc: _cs) bc->apply(grid, time); }
};

#endif //AEROHPC_A_BOUNDARIES_HPP
