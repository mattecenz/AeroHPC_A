#ifndef AEROHPC_A_BOUNDARYCONDITION_H
#define AEROHPC_A_BOUNDARYCONDITION_H

#include "Grid.hpp"
#include <functional>
#include <utility>

template<Addressing_T A>
class Condition {

public:

    /**
      * Typedef shortening lambda definition of mapping function
      */
    typedef std::function<void(Grid<A> &, const TFunction &, Real time)> Mapper;


    Condition() = delete;

    /**
     * Construct a boundary condition given the mapping function and the characteristic spatial function
     */
    Condition(const Mapper &mapper, TFunction function) : _mapper(mapper), _function(std::move(function)) {}

    /**
     * Apply the boundary condition onto the given grid
     */
    void apply(Grid<A> &grid, Real time) const { _mapper(grid, _function, time); }

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
