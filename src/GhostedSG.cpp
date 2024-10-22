#include "GhostedSG.hpp"

template<>
Real &GhostedSG<Addressing_T::STANDARD>::operator()(Component c, size_t i, size_t j, size_t k) {
    return _entries[c + ((i + _px) + ((j + _py) + (k + _pz) * nx) * ny) * Component::N_COMPONENTS];
}

template<>
const Real &GhostedSG<Addressing_T::STANDARD>::operator()(Component c, size_t i, size_t j, size_t k) const {
    return _entries[c + ((i + _px) + ((j + _py) + (k + _pz) * nx) * ny) * Component::N_COMPONENTS];
}