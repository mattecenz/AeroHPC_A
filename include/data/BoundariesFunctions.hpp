#ifndef BOUNDARIESFUNCTIONS_HPP
#define BOUNDARIESFUNCTIONS_HPP

#include "Traits.hpp"

#define BOUNDARY_DIRICHLET 0
#define BOUNDARY_NEUMANN 1

class BoundariesFunctions{
public:
    const int northType, southType, eastType, westType, frontType, backType;
    const TFunction *northFunction, *southFunction, *eastFunction, *westFunction, *frontFunction, *backFunction;

    BoundariesFunctions(
        const int northType, const TFunction *northFunction,
        const int southType, const TFunction *southFunction,
        const int eastType, const TFunction *eastFunction,
        const int westType, const TFunction *westFunction,
        const int frontType, const TFunction *frontFunction,
        const int backType, const TFunction *backFunction
    ) : northType(northType), southType(southType),
        eastType(eastType), westType(westType),
        frontType(frontType), backType(backType),
        northFunction(northFunction), southFunction(southFunction),
        eastFunction(eastFunction), westFunction(westFunction),
        frontFunction(frontFunction), backFunction(backFunction){};
};


#endif //BOUNDARIESFUNCTIONS_HPP
