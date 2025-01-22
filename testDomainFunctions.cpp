#include "Traits.hpp"
#include "data/DomainData.hpp"

TFunction exactU = [](const Real x, const Real y, const Real z, const Real t) -> Real{
  return std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
};

TFunction exactV = [](const Real x, const Real y, const Real z, const Real t) -> Real{
  return std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
};

TFunction exactW = [](const Real x, const Real y, const Real z, const Real t) -> Real{
  return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
};

TFunction exactP = [](const Real x, const Real y, const Real z, const Real t) -> Real{
  return std::cos(x)*std::cos(y)*std::cos(z)*std::sin(t);
};

SpatialFunctionsCollect northC(exactU, exactV, exactW, exactP);
SpatialFunctionsCollect southC(exactU, exactV, exactW, exactP);
SpatialFunctionsCollect backC(exactU, exactV, exactW, exactP);
SpatialFunctionsCollect frontC(exactU, exactV, exactW, exactP);
SpatialFunctionsCollect eastC(exactU, exactV, exactW, exactP);
SpatialFunctionsCollect westC(exactU, exactV, exactW, exactP);
SpatialFunctionsCollect domain(exactU, exactV, exactW, exactP);


DomainData testDomainData(BOUNDARY_DIRICHLET, northC,
                          BOUNDARY_DIRICHLET, southC,
                          BOUNDARY_DIRICHLET, eastC,
                          BOUNDARY_DIRICHLET, westC,
                          BOUNDARY_DIRICHLET, frontC,
                          BOUNDARY_DIRICHLET, backC,
                          domain);

/*
#include "Traits.hpp"
#include "data/DomainData.hpp"

TFunction zeroFunction = [](const Real x, const Real y, const Real z, const Real t) -> Real{
  return 0.0;
};

TFunction oneFunction = [](const Real x, const Real y, const Real z, const Real t) -> Real{
  return 1.0;
};

SpatialFunctionsCollect northC(zeroFunction, zeroFunction, zeroFunction, zeroFunction);
SpatialFunctionsCollect southC(zeroFunction, zeroFunction, zeroFunction, zeroFunction);
SpatialFunctionsCollect backC(zeroFunction, oneFunction, zeroFunction, zeroFunction);
SpatialFunctionsCollect frontC(zeroFunction, zeroFunction, zeroFunction, zeroFunction);
SpatialFunctionsCollect eastC(zeroFunction, zeroFunction, zeroFunction, zeroFunction);
SpatialFunctionsCollect westC(zeroFunction, zeroFunction, zeroFunction, zeroFunction);
SpatialFunctionsCollect domain(zeroFunction, zeroFunction, zeroFunction, zeroFunction);


DomainData testDomainData(BOUNDARY_DIRICHLET, northC,
                          BOUNDARY_DIRICHLET, southC,
                          BOUNDARY_DIRICHLET, eastC,
                          BOUNDARY_DIRICHLET, westC,
                          BOUNDARY_DIRICHLET, frontC,
                          BOUNDARY_DIRICHLET, backC,
                          domain);
*/