#ifndef AEROHPC_A_L2NORM_CALCULATOR_H
#define AEROHPC_A_L2NORM_CALCULATOR_H

#include "Traits.hpp"
#include "data/SolverData.hpp"
#include "utils/macroUtils.hpp"

inline void computeL2Norm(const Real *data, const Real time, Real &uErr, Real &pErr) {
    Real sumU = 0.0;
    Real sumP = 0.0;

    // Loop through the entire grid
    ITERATE_DOMAIN_VELOCITY(i, j, k, true, SKIP_U)
        const Real exactU = domData.domainF.UF(x + params.dX, y + params.dY2, z + params.dZ2, time);
        const Real gridU = U(data, i, j, k);
        const Real diffU = gridU - exactU;
        sumU += (diffU * diffU);
    ITERATE_DOMAIN_END()

    ITERATE_DOMAIN_VELOCITY(i, j, k, true, SKIP_V)
        const Real exactV = domData.domainF.VF(x + params.dX2, y + params.dY, z + params.dZ2, time);
        const Real gridV = V(data, i, j, k);
        const Real diffV = gridV - exactV;
        sumU += (diffV * diffV);
    ITERATE_DOMAIN_END()

    ITERATE_DOMAIN_VELOCITY(i, j, k, true, SKIP_W)
        const Real exactW = domData.domainF.WF(x + params.dX2, y + params.dY2, z + params.dZ, time);
        const Real gridW = W(data, i, j, k);
        const Real diffW = gridW - exactW;
        sumU += (diffW * diffW);
    ITERATE_DOMAIN_END()

    ITERATE_DOMAIN_PRESSURE(i, j, k, true)
        const Real exactP = domData.domainF.PF(x + params.dX2, y + params.dY2, z + params.dZ2, time);
        const Real gridP = P(data, i, j, k);
        const Real diffP = gridP - exactP;
        sumP += (diffP * diffP);
    ITERATE_DOMAIN_END()

    uErr = sumU * (params.dX * params.dY * params.dZ);
    pErr = sumP * (params.dX * params.dY * params.dZ);
}

#endif // AEROHPC_A_L2NORM_CALCULATOR_H
