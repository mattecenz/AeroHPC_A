#include "Traits.hpp"
#include "data/DomainData.hpp"
#include "data/SolverData.hpp"

//********************************* NON PERIODIC DOMAIN INFO **********************************************//
TFunction exactU = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y) * std::sin(2 * M_PI * z) * t;
};

TFunction exactV = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return std::cos(2 * M_PI * x) * std::sin(2 * M_PI * y) * std::sin(2 * M_PI * z) * t;
};

TFunction exactW = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return 2 * std::cos(2 * M_PI * x) * std::cos(2 * M_PI * y) * std::cos(2 * M_PI * z) * t;
};

TFunction exactP = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return std::cos(2 * M_PI * x) * std::cos(2 * M_PI * y) * std::sin(2 * M_PI * z) * t;
};

#define pw2(a) ((a)*(a))
TFunction forcingU = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return -2 * M_PI * pow(t, 2) * sin(2 * M_PI * x) * pow(sin(2 * M_PI * y), 2) * pow(sin(2 * M_PI * z), 2) * cos(2 * M_PI * x)
           + 2 * M_PI * pow(t, 2) * sin(2 * M_PI * x) * pow(sin(2 * M_PI * z), 2) * cos(2 * M_PI * x) * pow(cos(2 * M_PI * y), 2)
           + 4 * M_PI * pow(t, 2) * sin(2 * M_PI * x) * cos(2 * M_PI * x) * pow(cos(2 * M_PI * y), 2) * pow(cos(2 * M_PI * z), 2)
           - 2 * M_PI * t * sin(2 * M_PI * x) * sin(2 * M_PI * z) * cos(2 * M_PI * y)
           + sin(2 * M_PI * x) * sin(2 * M_PI * z) * cos(2 * M_PI * y)
           + 12 * pow(M_PI, 2) * t * sin(2 * M_PI * x) * sin(2 * M_PI * z) * cos(2 * M_PI * y) / params.Re;
};

TFunction forcingV = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return -2 * M_PI * pow(t, 2) * pow(sin(2 * M_PI * x), 2) * sin(2 * M_PI * y) * pow(sin(2 * M_PI * z), 2) * cos(2 * M_PI * y)
           + 2 * M_PI * pow(t, 2) * sin(2 * M_PI * y) * pow(sin(2 * M_PI * z), 2) * pow(cos(2 * M_PI * x), 2) * cos(2 * M_PI * y)
           + 4 * M_PI * pow(t, 2) * sin(2 * M_PI * y) * pow(cos(2 * M_PI * x), 2) * cos(2 * M_PI * y) * pow(cos(2 * M_PI * z), 2)
           - 2 * M_PI * t * sin(2 * M_PI * y) * sin(2 * M_PI * z) * cos(2 * M_PI * x)
           + sin(2 * M_PI * y) * sin(2 * M_PI * z) * cos(2 * M_PI * x)
           + 12 * pow(M_PI, 2) * t * sin(2 * M_PI * y) * sin(2 * M_PI * z) * cos(2 * M_PI * x) / params.Re;
};

TFunction forcingW = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return -4 * M_PI * pow(t, 2) * pow(sin(2 * M_PI * x), 2) * sin(2 * M_PI * z) * pow(cos(2 * M_PI * y), 2) * cos(2 * M_PI * z)
           - 4 * M_PI * pow(t, 2) * pow(sin(2 * M_PI * y), 2) * sin(2 * M_PI * z) * pow(cos(2 * M_PI * x), 2) * cos(2 * M_PI * z)
           - 8 * M_PI * pow(t, 2) * sin(2 * M_PI * z) * pow(cos(2 * M_PI * x), 2) * pow(cos(2 * M_PI * y), 2) * cos(2 * M_PI * z)
           + 2 * M_PI * t * cos(2 * M_PI * x) * cos(2 * M_PI * y) * cos(2 * M_PI * z)
           + 2 * cos(2 * M_PI * x) * cos(2 * M_PI * y) * cos(2 * M_PI * z)
           + 24 * pow(M_PI, 2) * t * cos(2 * M_PI * x) * cos(2 * M_PI * y) * cos(2 * M_PI * z) / params.Re;
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
