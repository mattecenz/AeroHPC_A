#include "ForcingTerm.hpp"

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
    return -pow(sin(time), 2) * sin(x) * pow(sin(y), 2) * pow(sin(z), 2) * cos(x)
           + pow(sin(time), 2) * sin(x) * pow(sin(z), 2) * cos(x) * pow(cos(y), 2)
           + 2.0 * pow(sin(time), 2) * sin(x) * cos(x) * pow(cos(y), 2) * pow(cos(z), 2)
           + sin(x) * sin(z) * cos(time) * cos(y)
           - 3.0 * sin(time) * sin(x) * sin(z) * cos(y) / Re;
}

// Private method to compute Gy
Real ForcingTerm::computeGy(Real x, Real y, Real z) const {
    return -pow(sin(time), 2) * pow(sin(x), 2) * sin(y) * pow(sin(z), 2) * cos(y)
           + pow(sin(time), 2) * sin(y) * pow(sin(z), 2) * pow(cos(x), 2) * cos(y)
           + 2.0 * pow(sin(time), 2) * sin(y) * pow(cos(x), 2) * cos(y) * pow(cos(z), 2)
           + sin(y) * sin(z) * cos(time) * cos(x)
           - 3.0 * sin(time) * sin(y) * sin(z) * cos(x) / Re;
}

// Private method to compute Gz
Real ForcingTerm::computeGz(Real x, Real y, Real z) const {
    return -2.0 * pow(sin(time), 2) * pow(sin(x), 2) * sin(z) * pow(cos(y), 2) * cos(z)
           - 2.0 * pow(sin(time), 2) * pow(sin(y), 2) * sin(z) * pow(cos(x), 2) * cos(z)
           - 4.0 * pow(sin(time), 2) * sin(z) * pow(cos(x), 2) * pow(cos(y), 2) * cos(z)
           + 2.0 * cos(time) * cos(x) * cos(y) * cos(z)
           - 6.0 * sin(time) * cos(x) * cos(y) * cos(z) / Re;
}


