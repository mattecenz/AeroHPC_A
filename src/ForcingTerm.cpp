#include "ForcingTerm.hpp"

// Constructor implementation
ForcingTerm::ForcingTerm(int Re, double time){
    set_Re(Re);
    set_time(time);
}

//setters
void ForcingTerm::set_time(double currentTime){
    time = currentTime;
}

void ForcingTerm::set_Re(int ReNum){
    Re = ReNum;
}


//getters
int ForcingTerm::get_Re() const{
    return Re;
}

double ForcingTerm::get_time() const{
    return time;
}
  
// Method to compute the forcing term at a given point in space and time
std::vector<double> ForcingTerm::compute(double x, double y, double z) const {
    vector<double> g(3);

    g[0] = computeGx(x, y, z); //Gx
    g[1] = computeGy(x, y, z); //Gy
    g[2] = computeGz(x, y, z); //Gz

    return g;
}

// Private method to compute Gx
double ForcingTerm::computeGx(double x, double y, double z) const {
    return ( sin(x)*cos(y)*sin(z)*(cos(get_time()) + (3.0/get_Re())*sin(get_time()))) +
    + (sin(x)*cos(x)*pow(sin(get_time()),2))*(pow(cos(y), 2)*pow(sin(z), 2)
    - pow(sin(y), 2)*pow(sin(z), 2) + 2*pow(cos(y),2)*pow(cos(z),2));    
}

// Private method to compute Gy
double ForcingTerm::computeGy(double x, double y, double z) const {
    return ( cos(x)*sin(y)*sin(z)*(cos(get_time()) + (3.0/get_Re())*sin(get_time()))) +
    + (cos(y)*sin(y)*pow(sin(get_time()),2))*(-pow(sin(x), 2)*pow(sin(z), 2)
    + pow(cos(x), 2)*pow(sin(z), 2) + 2*pow(cos(x),2)*pow(cos(z),2));    
}

// Private method to compute Gz
double ForcingTerm::computeGz(double x, double y, double z) const {
    return ( 2*cos(x)*cos(y)*cos(z)*(cos(get_time()) + (3.0/get_Re())*sin(get_time()))) +
    - 2*sin(z)*cos(z)*sin(get_time())* (pow(sin(x), 2)*pow(cos(y), 2)
    + pow(cos(x), 2)*pow(sin(y), 2) + 2*pow(cos(x),2)*pow(cos(y),2));    
}


