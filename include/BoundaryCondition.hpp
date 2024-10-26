#ifndef AEROHPC_A_BOUNDARYCONDITION_H
#define AEROHPC_A_BOUNDARYCONDITION_H

#include "Traits.hpp"
#include "StaggeredGrid.hpp"
#include <functional>
#include <utility>

template<Addressing_T A>
class BoundaryCondition {

public:

    /**
      * Typedef shortening lambda definition of mapping function
      */
    typedef std::function<void(StaggeredGrid<A> &, const Function &, Real time)> Mapper;


    BoundaryCondition() = delete;

    /**
     * Construct a boundary condition given the mapping function and the characteristic spatial function
     */
    BoundaryCondition(const Mapper &mapper, Function function) : _mapper(mapper), _function(std::move(function)) {}

    /**
     * Apply the boundary condition onto the given grid
     */
    void apply(StaggeredGrid<A> &grid, Real time) const { _mapper(grid, _function, time); }

private:
    /**
     * The BC mapping function
     */
    const Mapper _mapper;

    /**
     * The BC characteristic spatial function
     */
    const Function _function;
};

#endif //AEROHPC_A_BOUNDARYCONDITION_H
