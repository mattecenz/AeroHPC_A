#ifndef BOUNDARIESFUNCTIONS_HPP
#define BOUNDARIESFUNCTIONS_HPP

#include "Traits.hpp"

#define BOUNDARY_DIRICHLET 0
#define BOUNDARY_NEUMANN 1

class BoundaryFunction {
public:
    const TFunction *UF, *VF, *WF, *PF;

    BoundaryFunction(const TFunction *UF, const TFunction *VF,
                     const TFunction *WF, const TFunction *PF)
        : UF(UF), VF(VF), WF(WF), PF(PF) {
    }
};

class BoundariesFunctions {
public:
    const int northType, southType, eastType, westType, frontType, backType;
    const BoundaryFunction &northF, &southF, &eastF, &westF, &frontF, &backF;

    BoundariesFunctions(
        const int northType, const BoundaryFunction &northF,
        const int southType, const BoundaryFunction &southF,
        const int eastType, const BoundaryFunction &eastF,
        const int westType, const BoundaryFunction &westF,
        const int frontType, const BoundaryFunction &frontF,
        const int backType, const BoundaryFunction &backF
    ) : northType(northType), southType(southType),
        eastType(eastType), westType(westType),
        frontType(frontType), backType(backType),
        northF(northF), southF(southF),
        eastF(eastF), westF(westF),
        frontF(frontF), backF(backF) {
    };
};


#endif //BOUNDARIESFUNCTIONS_HPP
