#ifndef RKDATA_HPP
#define RKDATA_HPP

#include "Traits.hpp"
#include "Parameters.hpp"

#define VELOCITY_COMPONENTS 3

class RKData {
public:
    Real *model_data = nullptr;

    Real *buffer_data = nullptr;

    Real *rhs_data = nullptr;

    explicit RKData(const Parameters &params) {
        model_data = new Real[(VELOCITY_COMPONENTS + 1)  * params.grid_gndim];

        buffer_data = new Real[(VELOCITY_COMPONENTS + 1) * params.grid_gndim];

        rhs_data = new Real[(VELOCITY_COMPONENTS + 1) * params.grid_ndim];
    };

#define destroy(data) delete [] data; data=nullptr;

    ~RKData() {
        destroy(model_data);
        destroy(buffer_data);
        destroy(rhs_data);
    }
};

#endif //RKDATA_HPP
