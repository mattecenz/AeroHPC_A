#ifndef BOUNDARIESFUNCTIONS_HPP
#define BOUNDARIESFUNCTIONS_HPP

#include "data/SpatialFunctionsCollect.hpp"

/**
 * This class is a collector of information about domain and boundaries
 */
class DomainData {
public:
    /// This flags define whether the condition on the boundary is a Dirichlet or a Neumann one
#define BOUNDARY_DIRICHLET 0
#define BOUNDARY_NEUMANN 1
    const int northType, southType, eastType, westType, frontType, backType;

    /// Collections of functions that defines exact values of U, V, W and P on boundaries
    /// These collections are used only if the corresponding boundary type is Dirichlet
    const SpatialFunctionsCollect &northBF, &southBF, &eastBF, &westBF, &frontBF, &backBF;

    /// Collection of functions that defines exact values of U, V, W and P on all domain
    /// This collection is used only to check L2Norm on domain (so only for test purposes)
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
