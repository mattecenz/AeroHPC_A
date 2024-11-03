#ifndef AEROHPC_A_BOUNDARYCONDITION_H
#define AEROHPC_A_BOUNDARYCONDITION_H

#include "Grid.hpp"
#include <utility>

class Condition {

public:

    /**
      * Typedef shortening lambda definition of mapping function
      */
    typedef void (*Mapper)(Grid &, const TFunction &, Real);


    Condition() = delete;

    /**
     * Construct a boundary condition given the mapping function and the characteristic spatial function
     */
    Condition(const Mapper &mapper, TFunction function) : _mapper(mapper), _function(std::move(function)) {}

    /**
     * Apply the boundary condition onto the given grid
     */
    void apply(Grid &grid, const Real time) const { _mapper(grid, _function, time); }

private:
    /**
     * The BC mapping function
     */
    const Mapper _mapper;

    /**
     * The BC characteristic spatial function
     */
    const TFunction _function;
};

#endif //AEROHPC_A_BOUNDARYCONDITION_H
