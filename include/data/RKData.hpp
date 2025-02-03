#ifndef RKDATA_HPP
#define RKDATA_HPP

#include "Traits.hpp"
#include "data/Parameters.hpp"

class RKData {
public:
    /// Actual model buffer
    Real *model_data = nullptr;

    /// Alternative model buffer for RK
    Real *buffer_data = nullptr;

    /// RHS buffer for RK (velocity components)
    /// Base buffer of PS (pressure component)
    Real *rhs_data = nullptr;

    explicit RKData(const Parameters &params) {
        // Model buffers need ghosts
        model_data = new Real[(VELOCITY_COMPONENTS + 1)  * params.grid_gndim];
        buffer_data = new Real[(VELOCITY_COMPONENTS + 1) * params.grid_gndim];

        // RHS and Base Buffer don't need ghosts
        rhs_data = new Real[(VELOCITY_COMPONENTS + 1) * params.grid_ndim];
    }

    ~RKData() {
#define destroyData(data) delete [] data; data=nullptr;
        destroyData(model_data);
        destroyData(buffer_data);
        destroyData(rhs_data);
#undef destroyData
    }
};

#endif //RKDATA_HPP
