#include "ForcingTerm.hpp"

// Constructor implementation
ForcingTerm::ForcingTerm(Real Re, double time){
    set_Re(Re);
    set_time(time);
}

//setters
void ForcingTerm::set_time(double currentTime){
    time = currentTime;
}

void ForcingTerm::set_Re(Real ReNum){
    Re = ReNum;
}

void ForcingTerm::update_time(double extraTime){
    time += extraTime;
}

//getters
Real ForcingTerm::get_Re() const{
    return Re;
}


double ForcingTerm::get_time() const{
    return time;
}
  
// Method to compute the forcing term at a given point in space and time
double ForcingTerm::compute(double x, double y, double z, Component c) const{
    
    switch (c)
    {
    case Component::U:
            return computeGx(x, y, z, get_time(), get_Re()); //Gx
            break;
    case Component::V:
            return computeGy(x, y, z, get_time(), get_Re()); //Gy
            break;
    case Component::W:
            return computeGz(x, y, z, get_time(), get_Re()); //Gz
            break;
    }
    return 0.0;
}

// Private method to compute Gx
double ForcingTerm::computeGx(double x, double y, double z, double t, Real Re) const {
    return -pow(sin(t), 2) * sin(x) * pow(sin(y), 2) * pow(sin(z), 2) * cos(x) 
           + pow(sin(t), 2) * sin(x) * pow(sin(z), 2) * cos(x) * pow(cos(y), 2) 
           + 2 * pow(sin(t), 2) * sin(x) * cos(x) * pow(cos(y), 2) * pow(cos(z), 2) 
           + sin(x) * sin(z) * cos(t) * cos(y) 
           - 3 * sin(t) * sin(x) * sin(z) * cos(y) / Re;  
}

// Private method to compute Gy
double ForcingTerm::computeGy(double x, double y, double z, double t, Real Re) const {
    return -pow(sin(t), 2) * pow(sin(x), 2) * sin(y) * pow(sin(z), 2) * cos(y) 
           + pow(sin(t), 2) * sin(y) * pow(sin(z), 2) * pow(cos(x), 2) * cos(y) 
           + 2 * pow(sin(t), 2) * sin(y) * pow(cos(x), 2) * cos(y) * pow(cos(z), 2) 
           + sin(y) * sin(z) * cos(t) * cos(x) 
           - 3 * sin(t) * sin(y) * sin(z) * cos(x) / Re;
}

// Private method to compute Gz
double ForcingTerm::computeGz(double x, double y, double z, double t, Real Re) const {
    return -2 * pow(sin(t), 2) * pow(sin(x), 2) * sin(z) * pow(cos(y), 2) * cos(z) 
           - 2 * pow(sin(t), 2) * pow(sin(y), 2) * sin(z) * pow(cos(x), 2) * cos(z) 
           - 4 * pow(sin(t), 2) * sin(z) * pow(cos(x), 2) * pow(cos(y), 2) * cos(z) 
           + 2 * cos(t) * cos(x) * cos(y) * cos(z) 
           - 6 * sin(t) * cos(x) * cos(y) * cos(z) / Re;  
}


