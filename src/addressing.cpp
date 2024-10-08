#ifndef ADDRESSING_HPP
#define ADDRESSING_HPP

#include <vector>
#include <cstdio>

/**
 * Enum for easy modification of the class addressing
 */
enum class ADDRESSING{
    STANDARD=0
};

/**
 * 
 * Class for handling a 3D Matrix addressing.
 * 
 * Used as layer of abstraction in order to separate the addressing complexity
 * from the maths developement.
 * 
 * For the moment the ordering is: first x-y plane as if it were a rowmajor matrix then move on the z axis
 * 
 */
template <typename Real, ADDRESSING Order=ADDRESSING::STANDARD>
class Tensor{

public:

    // Constructors

    /**
     * Constructor of a nx * ny * nz tensor with initial value equal to zero
     */
    Tensor(const size_t nx, const size_t ny, const size_t nz):
        _nx(nx),
        _ny(ny),
        _nz(nz)
        {
            auto dim=nx*ny*nz;
            _entries.reserve(dim);
            for(size_t i=0;i<dim;++i){
                _entries.emplace_back(0.);
            }
        }

    /**
     * Direct accessing operator to an entry
     */
    Real& operator[](const size_t index){
        return _entries[index];
    }

    /**
     * Direct accessing operator for read-only operations
     */
    const Real& operator[](const size_t index) const{
        return _entries[index];
    }

    /**
     * Operator that accesses the memory using a 3D view of the object
     */
    Real& operator()(const size_t i, const size_t j, const size_t k){
        if constexpr(Order=ADDRESSING::STANDARD)
            return i + j*_ny + k*_nx*_ny;
    }

    /**
     * Operator that accesses 3D view for read-only operations
     */
    const Real& operator()(const size_t i, const size_t j, const size_t k) const{
        if constexpr(Order=ADDRESSING::STANDARD)
            return i + j*_nx + k*_nx*_ny;
    }

private:

    std::vector<Real> _entries;

    size_t _nx;
    size_t _ny;
    size_t _nz;

};

#endif