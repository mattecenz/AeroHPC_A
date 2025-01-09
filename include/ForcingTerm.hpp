#ifndef FORCING_TERM_HPP
#define FORCING_TERM_HPP

#include <array>
#include <cmath>
#include "Traits.hpp"

class ForcingTerm {

public:

    ForcingTerm(Real Re, Real time) : Re(Re), time(time){}

    void set_time(Real currentTime){
        time = currentTime;
    }

    void update_time(Real extraTime){
        time += extraTime;
    }

    Real get_time() const { return time; };

    Real computeGx(Real x, Real y, Real z) const;
    Real computeGy(Real x, Real y, Real z) const;
    Real computeGz(Real x, Real y, Real z) const;

private:
    Real time;
    const Real Re;
};

#endif // FORCING_TERM_HPP
