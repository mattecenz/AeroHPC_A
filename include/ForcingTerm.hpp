#ifndef FORCING_TERM_HPP
#define FORCING_TERM_HPP

#include <array>
#include <cmath>
#include "Traits.hpp"

class ForcingTerm {
    const TFunction FU, FV, FW;

public:
    ForcingTerm(const TFunction &forcingU, const TFunction &forcingV, const TFunction &forcingW)
        : FU(forcingU), FV(forcingV), FW(forcingW) {
    }

    void set_time(const Real currentTime) {
        time = currentTime;
    }

    void update_time(const Real extraTime) {
        time += extraTime;
    }

    Real get_time() const { return time; };

    Real computeU(const Real x, const Real y, const Real z) const { return FU(x, y, z, time); }
    Real computeV(const Real x, const Real y, const Real z) const { return FV(x, y, z, time); }
    Real computeW(const Real x, const Real y, const Real z) const { return FW(x, y, z, time); }

private:
    Real time = 0.0;
};

#endif // FORCING_TERM_HPP
