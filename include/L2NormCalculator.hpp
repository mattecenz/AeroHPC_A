#ifndef AEROHPC_A_L2NORM_CALCULATOR_H
#define AEROHPC_A_L2NORM_CALCULATOR_H

#include "Traits.hpp"
#include "data/SolverData.hpp"
#include "utils/macroUtils.hpp"

inline void computeL2Norm(const Real *data, const Real time, Real &uErr, Real &pErr) {
    Real sumU = 0.0;
    Real sumP = 0.0;

    // Loop through the entire grid
    ITERATE_DOMAIN_VELOCITY(i, j, k, true, NO_SKIP)
        const Real exactU = domData.domainF.UF(x + params.dX, y + params.dY2, z + params.dZ2, time);
        const Real exactV = domData.domainF.VF(x + params.dX2, y + params.dY, z + params.dZ2, time);
        const Real exactW = domData.domainF.WF(x + params.dX2, y + params.dY2, z + params.dZ, time);

        const Real exactP = domData.domainF.PF(x + params.dX2, y + params.dY2, z + params.dZ2, time);

        // Access the computed grid components
        const Real gridU = U(data, i, j, k);
        const Real gridV = V(data, i, j, k);
        const Real gridW = W(data, i, j, k);

        const Real gridP = P(data, i, j, k);

        // Calculate the differences
        const Real diffU = gridU - exactU;
        const Real diffV = gridV - exactV;
        const Real diffW = gridW - exactW;

        const Real diffP = gridP - exactP;

        // Add the squares of the differences to sum
        sumU += (diffU * diffU) + (diffV * diffV) + (diffW * diffW);

        sumP += (diffP * diffP);
    ITERATE_DOMAIN_END()

    uErr = sumU * (params.dX * params.dY * params.dZ);
    pErr = sumP * (params.dX * params.dY * params.dZ);
}

#endif // AEROHPC_A_L2NORM_CALCULATOR_H
