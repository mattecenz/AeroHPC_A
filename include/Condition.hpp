#ifndef AEROHPC_A_BOUNDARYCONDITION_H
#define AEROHPC_A_BOUNDARYCONDITION_H

#include "GridData.hpp"
#include <utility>

class Condition {

public:
    /**
     * Apply the boundary condition onto the given grid at given time
     */
    virtual void apply(GridData &grid, Real time) = 0;
};

class PhysicalCondition : public Condition {

public:

    /**
      * Typedef shortening lambda definition of mapping function
      */
    typedef void (*Mapper)(GridData &grid, const Real time, const std::vector<TFunction> &functions);


    PhysicalCondition() = delete;

    /**
     * Construct a boundary condition given the mapping function and the characteristic spatial function
     */
    PhysicalCondition(const Mapper &mapper, const std::vector<TFunction> &functions) : _mapper(mapper), _functions(functions) {}

    /**
     * Apply the boundary condition onto the given grid
     */
    void apply(GridData &grid, const Real time) override { _mapper(grid, time, _functions); }

private:
    /**
     * The BC mapping function
     */
    const Mapper _mapper;

    /**
     * The BC characteristic spatial functions
     */
    const std::vector<TFunction> _functions;
};

#endif //AEROHPC_A_BOUNDARYCONDITION_H
