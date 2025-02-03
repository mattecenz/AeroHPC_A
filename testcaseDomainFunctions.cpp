#include "Traits.hpp"
#include "data/DomainData.hpp"

TFunction zeroFunction = [](const Real x, const Real y, const Real z, const Real t) -> Real{
    return 0.0;
};

TFunction oneFunction = [](const Real x, const Real y, const Real z, const Real t) -> Real{
    return 1.0;
};

DomainData testcase1DomainData(DIRICHLET,  {zeroFunction, zeroFunction, zeroFunction, zeroFunction},
                               DIRICHLET, {zeroFunction, zeroFunction, zeroFunction, zeroFunction},
                               DIRICHLET, {zeroFunction, zeroFunction, zeroFunction, zeroFunction},
                               DIRICHLET, {zeroFunction, zeroFunction, zeroFunction, zeroFunction},
                               DIRICHLET, {zeroFunction, oneFunction, zeroFunction, zeroFunction},
                               DIRICHLET, {zeroFunction, zeroFunction, zeroFunction, zeroFunction},
                               {zeroFunction, zeroFunction, zeroFunction, zeroFunction});