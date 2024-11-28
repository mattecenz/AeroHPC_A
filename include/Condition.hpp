#ifndef AEROHPC_A_BOUNDARYCONDITION_H
#define AEROHPC_A_BOUNDARYCONDITION_H

#include "Grid.hpp"
#include <utility>

class Condition {

public:

    /**
      * Typedef shortening lambda definition of mapping function
      */
    typedef void (*Mapper)(Grid & grid, const Real time, const std::vector<TFunction>& functions);


    Condition() = delete;

    /**
     * Construct a boundary condition given the mapping function and the characteristic spatial function
     */
    Condition(const Mapper &mapper, const std::vector<TFunction>& functions) : _mapper(mapper), _functions(functions) {}

    /**
     * Apply the boundary condition onto the given grid
     */
    void apply(Grid &grid, const Real time) const { _mapper(grid,  time, _functions); }

private:
    /**
     * The BC mapping function
     */
    const Mapper _mapper;

    /**
     * The BC characteristic spatial function
     */
    const std::vector<TFunction> _functions;
};

#endif //AEROHPC_A_BOUNDARYCONDITION_H
