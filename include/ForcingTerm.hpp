#ifndef FORCING_TERM_HPP
#define FORCING_TERM_HPP

#include <array>
#include <cmath>  

using namespace std;

class ForcingTerm {

public:

    // Constructor 
    ForcingTerm(int Re, double time = 0.0);

    //setters
    void set_time(double time);
    void set_Re(int Re);

    //getters
    int get_Re() const;
    double get_time() const;

    // Method to compute the forcing term at a given point in space and time
    array<double,3> compute(double x, double y, double z) const;

    // Hfunctions to compute the forcing term in each direction.
    double computeGx(double x, double y, double z) const;
    double computeGy(double x, double y, double z) const;
    double computeGz(double x, double y, double z) const;

private:
    double time;
    int Re;
};

#endif // FORCING_TERM_HPP




