#ifndef FORCING_TERM_HPP
#define FORCING_TERM_HPP

#include <array>
#include <cmath>
#include "Grid.hpp"

using namespace std;

class ForcingTerm {

public:

    // Constructor 
    ForcingTerm(Real Re, Real time);

    //setters
    void set_time(Real time);

    //getters
    Real get_time() const;

    //Time update
    void update_time(Real extraTime);

    // Method to compute the forcing term at a given point in space and time
    Vector compute(Real x, Real y, Real z) const;

    // Hfunctions to compute the forcing term in each direction.
    Real computeGx(Real x, Real y, Real z) const;
    Real computeGy(Real x, Real y, Real z) const;
    Real computeGz(Real x, Real y, Real z) const;

private:
    Real time;
    Real Re;
};

#endif // FORCING_TERM_HPP
