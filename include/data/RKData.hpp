#ifndef RKDATA_HPP
#define RKDATA_HPP

#include "Traits.hpp"
#include "data/Parameters.hpp"

class RKData {
public:
    Real *model_data = nullptr;

    Real *buffer_data = nullptr;

    Real *rhs_data = nullptr;

    explicit RKData(const Parameters &params) {
        model_data = new Real[(VELOCITY_COMPONENTS + 1)  * params.grid_gndim];

        buffer_data = new Real[(VELOCITY_COMPONENTS + 1) * params.grid_gndim];

        rhs_data = new Real[(VELOCITY_COMPONENTS + 1) * params.grid_ndim];
    }

#define my_destroy(data) delete [] data; data=nullptr;

    ~RKData() {
        my_destroy(model_data);
        my_destroy(buffer_data);
        my_destroy(rhs_data);
    }
};

#endif //RKDATA_HPP
