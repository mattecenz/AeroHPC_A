#ifndef SPATIALFUNCTIONSCOLLECT_HPP
#define SPATIALFUNCTIONSCOLLECT_HPP

#include "Traits.hpp"

/**
 * This class contains the functions that defines exact values of Velocity and Pressure
 */
class SpatialFunctionsCollect {
public:
    const TFunction &UF, &VF, &WF, &PF;

    SpatialFunctionsCollect(const TFunction &UF, const TFunction &VF,
                     const TFunction &WF, const TFunction &PF)
        : UF(UF), VF(VF), WF(WF), PF(PF) {
    }
};

#endif //SPATIALFUNCTIONSCOLLECT_HPP
