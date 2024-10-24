#include "GhostedSG.hpp"

template<>
Real &GhostedSG<STANDARD>::operator()(Component c, size_t i, size_t j, size_t k) {
    return _entries[c + ((i + _px) + ((j + _py) + (k + _pz) * nx) * ny) * N_COMPONENTS];
}

template<>
const Real &GhostedSG<STANDARD>::operator()(Component c, size_t i, size_t j, size_t k) const {
    return _entries[c + ((i + _px) + ((j + _py) + (k + _pz) * nx) * ny) * N_COMPONENTS];
}