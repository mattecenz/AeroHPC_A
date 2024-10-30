#include "ForcingTerm.hpp"

#define sqrp(val) val * val

// Constructor implementation
ForcingTerm::ForcingTerm(Real Re, Real time) : Re(Re), time(time){}

//setters
void ForcingTerm::set_time(Real currentTime){
    time = currentTime;
}

void ForcingTerm::update_time(Real extraTime){
    time += extraTime;
}

// TODO not used, since it's not optimized for staggering single components
// Method to compute the forcing term at a given point in space and time
Vector ForcingTerm::compute(Real x, Real y, Real z) const {
    return {
            computeGx(x,y,z),
            computeGy(x,y,z),
            computeGz(x,y,z)
    };
}

// Private method to compute Gx
Real ForcingTerm::computeGx(Real x, Real y, Real z) const {
    return -sqrp(sin(time)) * sin(x) * sqrp(sin(y)) * sqrp(sin(z)) * cos(x)
           + sqrp(sin(time)) * sin(x) * sqrp(sin(z)) * cos(x) * sqrp(cos(y))
           + 2.0 * sqrp(sin(time)) * sin(x) * cos(x) * sqrp(cos(y)) * sqrp(cos(z))
           + sin(x) * sin(z) * cos(time) * cos(y)
           - 3.0 * sin(time) * sin(x) * sin(z) * cos(y) / Re;
}

// Private method to compute Gy
Real ForcingTerm::computeGy(Real x, Real y, Real z) const {
    return -sqrp(sin(time)) * sqrp(sin(x)) * sin(y) * sqrp(sin(z)) * cos(y)
           + sqrp(sin(time)) * sin(y) * sqrp(sin(z)) * sqrp(cos(x)) * cos(y)
           + 2.0 * sqrp(sin(time)) * sin(y) * sqrp(cos(x)) * cos(y) * sqrp(cos(z))
           + sin(y) * sin(z) * cos(time) * cos(x)
           - 3.0 * sin(time) * sin(y) * sin(z) * cos(x) / Re;
}

// Private method to compute Gz
Real ForcingTerm::computeGz(Real x, Real y, Real z) const {
    return -2.0 * sqrp(sin(time)) * sqrp(sin(x)) * sin(z) * sqrp(cos(y)) * cos(z)
           - 2.0 * sqrp(sin(time)) * sqrp(sin(y)) * sin(z) * sqrp(cos(x)) * cos(z)
           - 4.0 * sqrp(sin(time)) * sin(z) * sqrp(cos(x)) * sqrp(cos(y)) * cos(z)
           + 2.0 * cos(time) * cos(x) * cos(y) * cos(z)
           - 6.0 * sin(time) * cos(x) * cos(y) * cos(z) / Re;
}


