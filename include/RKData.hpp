#ifndef RKDATA_HPP
#define RKDATA_HPP

#include "Traits.hpp"
#include "Parameters.hpp"

#define VELOCITY_COMPONENTS 3

class RKData {
public:
    Real *velocity_data = nullptr;
    Real *pressure_data = nullptr;

    Real *velocity_buffer_data = nullptr;
    Real *pressure_buffer_data = nullptr;

    Real *rhs_buffer_data = nullptr;

    explicit RKData(const Parameters &params) {
        velocity_data = new Real[VELOCITY_COMPONENTS * params.grid_gndim];
        pressure_data = new Real[params.grid_gndim];

        velocity_buffer_data = new Real[VELOCITY_COMPONENTS * params.grid_gndim];
        pressure_buffer_data = new Real[params.grid_ndim];

        rhs_buffer_data = new Real[VELOCITY_COMPONENTS * params.grid_ndim];
    };

#define destroy(data) delete [] data; data=nullptr;

    ~RKData() {
        destroy(velocity_data);
        destroy(pressure_data);

        destroy(velocity_buffer_data);
        destroy(pressure_buffer_data);

        destroy(rhs_buffer_data);
    }
};

#endif //RKDATA_HPP
