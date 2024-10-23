#ifndef FORCING_TERM_HPP
#define FORCING_TERM_HPP

#include <array>
#include <cmath>  
#include <StaggeredGrid.hpp>

using namespace std;

class ForcingTerm {

public:

    // Constructor 
    ForcingTerm(Real Re, double time = 0.0);

    //setters
    void set_time(double time);
    void set_Re(Real Re);

    //getters
    Real get_Re() const;
    double get_time() const;

    //Time update
    void ForcingTerm::update_time(double extraTime);

    // Method to compute the forcing term at a given point in space and time
    double compute(double x, double y, double z, Component c) const;

    // Hfunctions to compute the forcing term in each direction.
    double computeGx(double x, double y, double z, double t, Real Re) const;
    double computeGy(double x, double y, double z, double t, Real Re) const;
    double computeGz(double x, double y, double z, double t, Real Re) const;

private:
    double time;
    Real Re;
};

#endif // FORCING_TERM_HPP
