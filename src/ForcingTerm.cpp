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
    return (Re * (-2 * pow(sin(time), 2) * pow(sin(y), 2) * cos(x) 
                  - pow(sin(time), 2) * pow(sin(z), 2) * cos(x) 
                  + 2 * pow(sin(time), 2) * cos(x) 
                  + sin(z) * cos(time) * cos(y)) 
            + 3 * sin(time) * sin(z) * cos(y)) * sin(x) / Re;
}

// Private method to compute Gy
Real ForcingTerm::computeGy(Real x, Real y, Real z) const {
    return (Re * (-2 * pow(sin(time), 2) * pow(sin(x), 2) * cos(y) 
                  - pow(sin(time), 2) * pow(sin(z), 2) * cos(y) 
                  + 2 * pow(sin(time), 2) * cos(y) 
                  + sin(z) * cos(time) * cos(x)) 
            + 3 * sin(time) * sin(z) * cos(x)) * sin(y) / Re;
}

// Private method to compute Gz
Real ForcingTerm::computeGz(Real x, Real y, Real z) const {
    return 2 * (Re * (pow(sin(time), 2) * pow(sin(x), 2) * sin(z) 
                      + pow(sin(time), 2) * pow(sin(y), 2) * sin(z) 
                      - 2 * pow(sin(time), 2) * sin(z) 
                      + cos(time) * cos(x) * cos(y)) 
                + 3 * sin(time) * cos(x) * cos(y)) * cos(z) / Re;
}


