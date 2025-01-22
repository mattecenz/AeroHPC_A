#include "ForcingTerm.hpp"

#define pw2(a) ((a)*(a))

#if DISABLE_PRESSURE
Real ForcingTerm::computeGx(Real x, Real y, Real z) const {
    return - pw2(sin(time)) * sin(x) * pw2(sin(y)) * pw2(sin(z)) * cos(x)
           + pw2(sin(time)) * sin(x) * pw2(sin(z)) * cos(x) * pw2(cos(y))
           + 2.0 * pw2(sin(time)) * sin(x) * cos(x) * pw2(cos(y)) * pw2(cos(z))
           + sin(x) * sin(z) * cos(time) * cos(y)
           + 3.0 * sin(time) * sin(x) * sin(z) * cos(y) / Re;
}

Real ForcingTerm::computeGy(Real x, Real y, Real z) const {
    return - pw2(sin(time)) * pw2(sin(x)) * sin(y) * pw2(sin(z)) * cos(y)
           + pw2(sin(time)) * sin(y) * pw2(sin(z)) * pw2(cos(x)) * cos(y)
           + 2.0 * pw2(sin(time)) * sin(y) * pw2(cos(x)) * cos(y) * pw2(cos(z))
           + sin(y) * sin(z) * cos(time) * cos(x)
           + 3.0 * sin(time) * sin(y) * sin(z) * cos(x) / Re;
}

Real ForcingTerm::computeGz(Real x, Real y, Real z) const {
    return - 2.0 * pw2(sin(time)) * pw2(sin(x)) * sin(z) * pw2(cos(y)) * cos(z)
           - 2.0 * pw2(sin(time)) * pw2(sin(y)) * sin(z) * pw2(cos(x)) * cos(z)
           - 4.0 * pw2(sin(time)) * sin(z) * pw2(cos(x)) * pw2(cos(y)) * cos(z)
           + 2.0 * cos(time) * cos(x) * cos(y) * cos(z)
           + 6.0 * sin(time) * cos(x) * cos(y) * cos(z) / Re;
}

#else
Real ForcingTerm::computeGx(Real x, Real y, Real z) const {
    return -2.0 * M_PI * sin(time) * sin(2.0 * M_PI * x) * cos(2.0 * M_PI * y) * cos(2.0 * M_PI * z)
           - pw2(sin(time)) * sin(x) * pw2(sin(y)) * pw2(sin(z)) * cos(x)
           + pw2(sin(time)) * sin(x) * pw2(sin(z)) * cos(x) * pw2(cos(y))
           + 2.0 * pw2(sin(time)) * sin(x) * cos(x) * pw2(cos(y)) * pw2(cos(z))
           + sin(x) * sin(z) * cos(time) * cos(y)
           + 3.0 * sin(time) * sin(x) * sin(z) * cos(y) / Re;
//    return 0.0;
}

Real ForcingTerm::computeGy(Real x, Real y, Real z) const {
    return -2.0 * M_PI * sin(time) * sin(2.0 * M_PI * y) * cos(2.0 * M_PI * x) * cos(2.0 * M_PI * z)
           - pw2(sin(time)) * pw2(sin(x)) * sin(y) * pw2(sin(z)) * cos(y)
           + pw2(sin(time)) * sin(y) * pw2(sin(z)) * pw2(cos(x)) * cos(y)
           + 2.0 * pw2(sin(time)) * sin(y) * pw2(cos(x)) * cos(y) * pw2(cos(z))
           + sin(y) * sin(z) * cos(time) * cos(x)
           + 3.0 * sin(time) * sin(y) * sin(z) * cos(x) / Re;
//    return 0.0;
}

Real ForcingTerm::computeGz(Real x, Real y, Real z) const {
    return -2.0 * M_PI * sin(time) * sin(2.0 * M_PI * z) * cos(2.0 * M_PI * x) * cos(2.0 * M_PI * y)
           - 2.0 * pw2(sin(time)) * pw2(sin(x)) * sin(z) * pw2(cos(y)) * cos(z)
           - 2.0 * pw2(sin(time)) * pw2(sin(y)) * sin(z) * pw2(cos(x)) * cos(z)
           - 4.0 * pw2(sin(time)) * sin(z) * pw2(cos(x)) * pw2(cos(y)) * cos(z)
           + 2.0 * cos(time) * cos(x) * cos(y) * cos(z)
           + 6.0 * sin(time) * cos(x) * cos(y) * cos(z) / Re;
//    return 0.0;
}
#endif
