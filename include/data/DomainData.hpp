#ifndef BOUNDARIESFUNCTIONS_HPP
#define BOUNDARIESFUNCTIONS_HPP

#include "data/SpatialFunctionsCollect.hpp"

#define BOUNDARY_DIRICHLET 0
#define BOUNDARY_NEUMANN 1

class DomainData {
public:
    const int northType, southType, eastType, westType, frontType, backType;
    const SpatialFunctionsCollect &northBF, &southBF, &eastBF, &westBF, &frontBF, &backBF;
    const SpatialFunctionsCollect &domainF;

    DomainData() = delete;

    DomainData(
        const int northType, const SpatialFunctionsCollect &northBF,
        const int southType, const SpatialFunctionsCollect &southBF,
        const int eastType, const SpatialFunctionsCollect &eastBF,
        const int westType, const SpatialFunctionsCollect &westBF,
        const int frontType, const SpatialFunctionsCollect &frontBF,
        const int backType, const SpatialFunctionsCollect &backBF,
        const SpatialFunctionsCollect &domainF
    ) : northType(northType), southType(southType),
        eastType(eastType), westType(westType),
        frontType(frontType), backType(backType),
        northBF(northBF), southBF(southBF),
        eastBF(eastBF), westBF(westBF),
        frontBF(frontBF), backBF(backBF),
        domainF(domainF) {
    }
};


#endif //BOUNDARIESFUNCTIONS_HPP
