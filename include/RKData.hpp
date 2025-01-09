#ifndef RKDATA_HPP
#define RKDATA_HPP

#include "Traits.hpp"
#include "Parameters.hpp"

#define VELOCITY_COMPONENTS 3

class RKData{
public:
    Real *velocity_data = nullptr;
    Real *pressure_data = nullptr;

    Real *velocity_buffer = nullptr;
    Real *pressure_buffer = nullptr;

    Real *rhs_buffer = nullptr;

    explicit RKData(const Parameters &params){
        velocity_data = new Real[VELOCITY_COMPONENTS * params.loc_gnX * params.loc_gnY * params.loc_gnZ];

        pressure_data = new Real[params.loc_gnX * params.loc_gnY * params.loc_gnZ];

        velocity_buffer = new Real[VELOCITY_COMPONENTS * params.loc_gnX * params.loc_gnY * params.loc_gnZ];

        pressure_buffer = new Real[params.loc_gnX * params.loc_gnY * params.loc_gnZ];

        rhs_buffer = new Real[3 * params.loc_nX * params.loc_nY * params.loc_nZ];
    };

    ~RKData(){
        delete[] velocity_data;
        delete[] pressure_data;
        delete[] velocity_buffer;
        delete[] pressure_buffer;
        delete[] rhs_buffer;
    }
};

#endif //RKDATA_HPP
