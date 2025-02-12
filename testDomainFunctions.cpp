#include "Traits.hpp"
#include "data/DomainData.hpp"
#include "data/SolverData.hpp"

//********************************* NON PERIODIC DOMAIN INFO **********************************************//
TFunction exactU = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
};

TFunction exactV = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
};

TFunction exactW = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
};

TFunction exactP = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return std::cos(x) * std::cos(y) * std::sin(z) * std::sin(t);
};

#define pw2(a) ((a)*(a))
TFunction forcingU = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return -2.0 * M_PI * sin(t) * sin(2.0 * M_PI * x) * cos(2.0 * M_PI * y) * cos(2.0 * M_PI * z)
           - pw2(sin(t)) * sin(x) * pw2(sin(y)) * pw2(sin(z)) * cos(x)
           + pw2(sin(t)) * sin(x) * pw2(sin(z)) * cos(x) * pw2(cos(y))
           + 2.0 * pw2(sin(t)) * sin(x) * cos(x) * pw2(cos(y)) * pw2(cos(z))
           + sin(x) * sin(z) * cos(t) * cos(y)
           + 3.0 * sin(t) * sin(x) * sin(z) * cos(y) / params.Re;
};

TFunction forcingV = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return -2.0 * M_PI * sin(t) * sin(2.0 * M_PI * y) * cos(2.0 * M_PI * x) * cos(2.0 * M_PI * z)
           - pw2(sin(t)) * pw2(sin(x)) * sin(y) * pw2(sin(z)) * cos(y)
           + pw2(sin(t)) * sin(y) * pw2(sin(z)) * pw2(cos(x)) * cos(y)
           + 2.0 * pw2(sin(t)) * sin(y) * pw2(cos(x)) * cos(y) * pw2(cos(z))
           + sin(y) * sin(z) * cos(t) * cos(x)
           + 3.0 * sin(t) * sin(y) * sin(z) * cos(x) / params.Re;
};

TFunction forcingW = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return -2.0 * M_PI * sin(t) * sin(2.0 * M_PI * z) * cos(2.0 * M_PI * x) * cos(2.0 * M_PI * y)
           - 2.0 * pw2(sin(t)) * pw2(sin(x)) * sin(z) * pw2(cos(y)) * cos(z)
           - 2.0 * pw2(sin(t)) * pw2(sin(y)) * sin(z) * pw2(cos(x)) * cos(z)
           - 4.0 * pw2(sin(t)) * sin(z) * pw2(cos(x)) * pw2(cos(y)) * cos(z)
           + 2.0 * cos(t) * cos(x) * cos(y) * cos(z)
           + 6.0 * sin(t) * cos(x) * cos(y) * cos(z) / params.Re;
};
#undef pw2

DomainData testDomainData(DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          {exactU, exactV, exactW, exactP},
                          {forcingU, forcingV, forcingW});
//*********************************************************************************************************//
