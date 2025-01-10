#ifndef SPATIALFUNCTIONSCOLLECT_HPP
#define SPATIALFUNCTIONSCOLLECT_HPP

#include "Traits.hpp"

class SpatialFunctionsCollect {
public:
    const TFunction *UF, *VF, *WF, *PF;

    SpatialFunctionsCollect(const TFunction *UF, const TFunction *VF,
                     const TFunction *WF, const TFunction *PF)
        : UF(UF), VF(VF), WF(WF), PF(PF) {
    }
};

#endif //SPATIALFUNCTIONSCOLLECT_HPP
