#include "Traits.hpp"
#include "data/DomainData.hpp"
#include "data/SolverData.hpp"
#include <cmath>

//********************************* NON PERIODIC DOMAIN INFO **********************************************//
TFunction exactU = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return sin(x) * cos(y) * sin(z) * sin(t);
};

TFunction exactV = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return cos(x) * sin(y) * sin(z) * sin(t);
};

TFunction exactW = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return 2 * cos(x) * cos(y) * cos(z) * sin(t);
};

TFunction exactP = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return cos(x) * cos(y) * sin(z) * sin(t);
};

TFunction forcingU = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return -pow(sin(t), 2) * sin(x) * pow(sin(y), 2) * pow(sin(z), 2) * cos(x) 
    + pow(sin(t), 2) * sin(x) * pow(sin(z), 2) * cos(x) * pow(cos(y), 2) 
    + 2 * pow(sin(t), 2) * sin(x) * cos(x) * pow(cos(y), 2) * pow(cos(z), 2) 
    + sin(x) * sin(z) * cos(t) * cos(y) 
    - sin(x) * sin(z) * cos(y) * sin(t)
    + 3 * sin(t) * sin(x) * sin(z) * cos(y) / params.Re;
};

TFunction forcingV = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return -pow(sin(t), 2) * pow(sin(x), 2) * sin(y) * pow(sin(z), 2) * cos(y) 
    + pow(sin(t), 2) * sin(y) * pow(sin(z), 2) * pow(cos(x), 2) * cos(y) 
    + 2 * pow(sin(t), 2) * sin(y) * pow(cos(x), 2) * cos(y) * pow(cos(z), 2) 
    + sin(y) * sin(z) * cos(t) * cos(x) 
    - sin(y) * sin(z) * cos(x) * sin(t)
    + 3 * sin(t) * sin(y) * sin(z) * cos(x) / params.Re;
};

TFunction forcingW = [](const Real x, const Real y, const Real z, const Real t) -> Real {
    return -2 * pow(sin(t), 2) * pow(sin(x), 2) * sin(z) * pow(cos(y), 2) * cos(z) 
    - 2 * pow(sin(t), 2) * pow(sin(y), 2) * sin(z) * pow(cos(x), 2) * cos(z) 
    - 4 * pow(sin(t), 2) * sin(z) * pow(cos(x), 2) * pow(cos(y), 2) * cos(z) 
    + 2 * cos(t) * cos(x) * cos(y) * cos(z) 
    + cos(x) * cos(y) * cos(z) * sin(t)
    + 6 * sin(t) * cos(x) * cos(y) * cos(z) / params.Re;
};

DomainData testDomainData(DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          {exactU, exactV, exactW, exactP},
                          {forcingU, forcingV, forcingW});
//*********************************************************************************************************//
