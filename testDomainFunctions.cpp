#include "Traits.hpp"
#include "data/DomainData.hpp"

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
    return std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
};

DomainData testDomainData(DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          DIRICHLET, {exactU, exactV, exactW, exactP},
                          {exactU, exactV, exactW, exactP});
