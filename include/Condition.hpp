#ifndef AEROHPC_A_BOUNDARYCONDITION_H
#define AEROHPC_A_BOUNDARYCONDITION_H

#include "Grid.hpp"
#include <utility>

class Condition {

public:
    /**
     * Apply the boundary condition onto the given grid
     */
    virtual void apply(Grid &grid, Real time) const = 0;
};

#endif //AEROHPC_A_BOUNDARYCONDITION_H
