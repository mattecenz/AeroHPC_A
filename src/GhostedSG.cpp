#include "GhostedSG.hpp"

template<>
Real &GhostedSG<STANDARD>::operator()(Component c, index_t i, index_t j, index_t k) {
    return _entries[c + ((i + _px) + ((j + _py) + (k + _pz) * (ny + 2*_py)) * (nx + 2*_px)) * N_COMPONENTS];
}

template<>
const Real &GhostedSG<STANDARD>::operator()(Component c, index_t i, index_t j, index_t k) const {
    return _entries[c + ((i + _px) + ((j + _py) + (k + _pz) * (ny + 2*_py)) * (nx + 2*_px)) * N_COMPONENTS];
}