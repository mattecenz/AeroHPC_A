
#include "Traits.hpp"
#include "data/DomainData.hpp"

TFunction exactU = [](Real x, Real y, Real z, Real t) -> Real{
  return 0;
};

TFunction exactV = [](Real x, Real y, Real z, Real t) -> Real{
  return 0;
};

TFunction exactW = [](Real x, Real y, Real z, Real t) -> Real{
  return 0;
};

TFunction exactP = [](Real x, Real y, Real z, Real t) -> Real{
  return 0;
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

