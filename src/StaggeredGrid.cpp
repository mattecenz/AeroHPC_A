#include <StaggeredGrid.hpp>

/**
 * Define standard addressing methods for the StaggeredGrid class
 */
template<>
Real &StaggeredGrid<Addressing_T::STANDARD>::operator()(const Component c,
                                                              const size_t i,
                                                              const size_t j,
                                                              const size_t k) {
    return _entries[c + (i + (j + k * nx) * ny) * Component::N_COMPONENTS];
}

/**
 * Define standard addressing methods for the StaggeredGrid class
 */
template<>
const Real &StaggeredGrid<Addressing_T::STANDARD>::operator()(const Component c,
                                                                    const size_t i,
                                                                    const size_t j,
                                                                    const size_t k) const {
    return _entries[c + (i + (j + k * nx) * ny) * Component::N_COMPONENTS];
}
