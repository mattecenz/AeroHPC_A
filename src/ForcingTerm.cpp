#include "ForcingTerm.hpp"

Real ForcingTerm::computeGx(Real x, Real y, Real z) const {
    return (Re * (-2 * pow(sin(time), 2) * pow(sin(y), 2) * cos(x)
                  - pow(sin(time), 2) * pow(sin(z), 2) * cos(x)
                  + 2 * pow(sin(time), 2) * cos(x)
                  + sin(z) * cos(time) * cos(y))
            + 3 * sin(time) * sin(z) * cos(y)) * sin(x) / Re;
}

Real ForcingTerm::computeGy(Real x, Real y, Real z) const {
    return (Re * (-2 * pow(sin(time), 2) * pow(sin(x), 2) * cos(y)
                  - pow(sin(time), 2) * pow(sin(z), 2) * cos(y)
                  + 2 * pow(sin(time), 2) * cos(y)
                  + sin(z) * cos(time) * cos(x))
            + 3 * sin(time) * sin(z) * cos(x)) * sin(y) / Re;
}

Real ForcingTerm::computeGz(Real x, Real y, Real z) const {
    return 2 * (Re * (pow(sin(time), 2) * pow(sin(x), 2) * sin(z)
                      + pow(sin(time), 2) * pow(sin(y), 2) * sin(z)
                      - 2 * pow(sin(time), 2) * sin(z)
                      + cos(time) * cos(x) * cos(y))
                + 3 * sin(time) * cos(x) * cos(y)) * cos(z) / Re;
}


