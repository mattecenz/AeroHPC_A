#ifndef AEROHPC_A_CONDITIONS_HPP
#define AEROHPC_A_CONDITIONS_HPP

#include "Condition.hpp"

class PhysicalCondition: public Condition {

    public:

        /**
          * Typedef shortening lambda definition of mapping function
          */
        typedef void (*Mapper)(Grid & grid, const Real time, const std::vector<TFunction>& functions);


        PhysicalCondition() = delete;

        /**
         * Construct a boundary condition given the mapping function and the characteristic spatial function
         */
        PhysicalCondition(const Mapper &mapper, const std::vector<TFunction>& functions) : _mapper(mapper), _functions(functions) {}

        /**
         * Apply the boundary condition onto the given grid
         */
        void apply(Grid &grid, const Real time) const override { _mapper(grid,  time, _functions); }

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

class GhostCondition: public Condition {
public:

    /**
      * Typedef shortening lambda definition of mapping function
      */
    typedef void (*Mapper)(Grid & grid, const Real time, const std::vector<DataFunction>& functions);


    GhostCondition() = delete;

    /**
     * Construct a boundary condition given the mapping function and the characteristic data function
     */
    GhostCondition(const Mapper &mapper, const std::vector<DataFunction>& functions) : _mapper(mapper), _functions(functions) {}

    /**
     * Apply the boundary condition onto the given grid
     */
    void apply(Grid &grid, const Real time) const override { _mapper(grid,  time, _functions); }

private:
    /**
     * The BC mapping function
     */
    const Mapper _mapper;

    /**
     * The BC characteristic data functions
     */
    const std::vector<DataFunction> _functions;
};


#endif //AEROHPC_A_CONDITIONS_HPP
