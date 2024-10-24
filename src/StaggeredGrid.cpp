#include <StaggeredGrid.hpp>

/**
 * Define standard addressing methods for the StaggeredGrid class
 */
template<>
Real &StaggeredGrid<STANDARD>::operator()(const Component c,
                                                              const index_t i,
                                                              const index_t j,
                                                              const index_t k) {
    return _entries[c + (i + (j + k * nx) * ny) * N_COMPONENTS];
}

/**
 * Define standard addressing methods for the StaggeredGrid class
 */
template<>
const Real &StaggeredGrid<STANDARD>::operator()(const Component c,
                                                                    const index_t i,
                                                                    const index_t j,
                                                                    const index_t k) const {
    return _entries[c + (i + (j + k * nx) * ny) * N_COMPONENTS];
}
