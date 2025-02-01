#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "data/SolverData.hpp"

#include "utils/macroUtils.hpp"
#include "utils/mathUtils.hpp"

namespace mu = mathUtils;

inline void interpolateData(const Real *data, const Real time) {
    // INTERPOLATE INTO DOMAIN (WITHOUT BOUNDARIES)
    ITERATE_DOMAIN_PHYSICAL(i, j, k, false)
        PU(interpData_ptr, i, j, k) = mu::interp_U_onGrid(data, i, j, k);
        PV(interpData_ptr, i, j, k) = mu::interp_V_onGrid(data, i, j, k);
        PW(interpData_ptr, i, j, k) = mu::interp_W_onGrid(data, i, j, k);
        PP(interpData_ptr, i, j, k) = mu::interp_P_onGrid(data, i, j, k);
    ITERATE_DOMAIN_END()

    // INTERPOLATE NORTH BOUNDARY (WITHOUT EDGES)
    const index_t north_j = params.phy_nY - 1;
    if (params.isOnTop && !params.periodicY) {
        const Real Y = real(north_j + params.st_nY) * params.dY + params.originY;
        ITERATE_XZ_PHYSICAL_FACE(i, k, true)
            PU(interpData_ptr, i, north_j, k) = domData.northBF.UF(X, Y, Z, time);
            PV(interpData_ptr, i, north_j, k) = domData.northBF.VF(X, Y, Z, time);
            PW(interpData_ptr, i, north_j, k) = domData.northBF.WF(X, Y, Z, time);
            PP(interpData_ptr, i, north_j, k) = domData.northBF.PF(X, Y, Z, time);
        ITERATE_FACE_END()
    } else {
        ITERATE_XZ_PHYSICAL_FACE(i, k, false)
            PU(interpData_ptr, i, north_j, k) = mu::interp_U_onGrid(data, i, north_j, k);
            PV(interpData_ptr, i, north_j, k) = mu::interp_V_onGrid(data, i, north_j, k);
            PW(interpData_ptr, i, north_j, k) = mu::interp_W_onGrid(data, i, north_j, k);
            PP(interpData_ptr, i, north_j, k) = mu::interp_P_onGrid(data, i, north_j, k);
        ITERATE_FACE_END()
    }

    // INTERPOLATE SOUTH BOUNDARY (WITHOUT EDGES)
    constexpr index_t south_j = 0;
    if (params.isOnBottom && !params.periodicY) {
        const Real Y = real(south_j + params.st_nY) * params.dY + params.originY;
        ITERATE_XZ_PHYSICAL_FACE(i, k, true)
            PU(interpData_ptr, i, south_j, k) = domData.southBF.UF(X, Y, Z, time);
            PV(interpData_ptr, i, south_j, k) = domData.southBF.VF(X, Y, Z, time);
            PW(interpData_ptr, i, south_j, k) = domData.southBF.WF(X, Y, Z, time);
            PP(interpData_ptr, i, south_j, k) = domData.southBF.PF(X, Y, Z, time);
        ITERATE_FACE_END()
    } else {
        ITERATE_XZ_PHYSICAL_FACE(i, k, false)
            PU(interpData_ptr, i, south_j, k) = mu::interp_U_onGrid(data, i, south_j, k);
            PV(interpData_ptr, i, south_j, k) = mu::interp_V_onGrid(data, i, south_j, k);
            PW(interpData_ptr, i, south_j, k) = mu::interp_W_onGrid(data, i, south_j, k);
            PP(interpData_ptr, i, south_j, k) = mu::interp_P_onGrid(data, i, south_j, k);
        ITERATE_FACE_END()
    }

    // INTERPOLATE EAST BOUNDARY (WITHOUT EDGES)
    const index_t east_k = params.phy_nZ - 1;
    if (params.isOnRight && !params.periodicZ) {
        const Real Z = real(east_k + params.st_nZ) * params.dZ + params.originZ;
        ITERATE_XY_PHYSICAL_FACE(i, j, true)
            PU(interpData_ptr, i, j, east_k) = domData.eastBF.UF(X, Y, Z, time);
            PV(interpData_ptr, i, j, east_k) = domData.eastBF.VF(X, Y, Z, time);
            PW(interpData_ptr, i, j, east_k) = domData.eastBF.WF(X, Y, Z, time);
            PP(interpData_ptr, i, j, east_k) = domData.eastBF.PF(X, Y, Z, time);
        ITERATE_FACE_END()
    } else {
        ITERATE_XY_PHYSICAL_FACE(i, j, false)
            PU(interpData_ptr, i, j, east_k) = mu::interp_U_onGrid(data, i, j, east_k);
            PV(interpData_ptr, i, j, east_k) = mu::interp_V_onGrid(data, i, j, east_k);
            PW(interpData_ptr, i, j, east_k) = mu::interp_W_onGrid(data, i, j, east_k);
            PP(interpData_ptr, i, j, east_k) = mu::interp_P_onGrid(data, i, j, east_k);
        ITERATE_FACE_END()
    }

    // INTERPOLATE WEST BOUNDARY (WITHOUT EDGES)
    constexpr index_t west_k = 0;
    if (params.isOnLeft && !params.periodicZ) {
        const Real Z = real(west_k + params.st_nZ) * params.dZ + params.originZ;
        ITERATE_XY_PHYSICAL_FACE(i, j, true)
            PU(interpData_ptr, i, j, west_k) = domData.westBF.UF(X, Y, Z, time);
            PV(interpData_ptr, i, j, west_k) = domData.westBF.VF(X, Y, Z, time);
            PW(interpData_ptr, i, j, west_k) = domData.westBF.WF(X, Y, Z, time);
            PP(interpData_ptr, i, j, west_k) = domData.westBF.PF(X, Y, Z, time);
        ITERATE_FACE_END()
    } else {
        ITERATE_XY_PHYSICAL_FACE(i, j, false)
            PU(interpData_ptr, i, j, west_k) = mu::interp_U_onGrid(data, i, j, west_k);
            PV(interpData_ptr, i, j, west_k) = mu::interp_V_onGrid(data, i, j, west_k);
            PW(interpData_ptr, i, j, west_k) = mu::interp_W_onGrid(data, i, j, west_k);
            PP(interpData_ptr, i, j, west_k) = mu::interp_P_onGrid(data, i, j, west_k);
        ITERATE_FACE_END()
    }

    // INTERPOLATE BACK BOUNDARY
    const index_t back_i = params.phy_nX - 1;
    if (!params.periodicX) {
        const Real X = real(back_i + params.st_nX) * params.dX + params.originX;
        for (index_t k = 0; k < params.phy_nZ; k++) {
            Real Z = real(k + params.st_nZ) * params.dZ + params.originZ;
            for (index_t j = 0; j < params.phy_nY; j++) {
                const Real Y = real(j + params.st_nY) * params.dY + params.originY;
                PU(interpData_ptr, back_i, j, k) = domData.westBF.UF(X, Y, Z, time);
                PV(interpData_ptr, back_i, j, k) = domData.westBF.VF(X, Y, Z, time);
                PW(interpData_ptr, back_i, j, k) = domData.westBF.WF(X, Y, Z, time);
                PP(interpData_ptr, back_i, j, k) = domData.westBF.PF(X, Y, Z, time);
            }
        }
    } else {
        for (index_t k = 0; k < params.phy_nZ; k++) {
            for (index_t j = 0; j < params.phy_nY; j++) {
                PU(interpData_ptr, back_i, j, k) = mu::interp_U_onGrid(data, back_i, j, k);
                PV(interpData_ptr, back_i, j, k) = mu::interp_V_onGrid(data, back_i, j, k);
                PW(interpData_ptr, back_i, j, k) = mu::interp_W_onGrid(data, back_i, j, k);
                PP(interpData_ptr, back_i, j, k) = mu::interp_P_onGrid(data, back_i, j, k);
            }
        }
    }

    // INTERPOLATE FRONT BOUNDARY
    constexpr index_t front_i = 0;
    if (!params.periodicX) {
        const Real X = real(front_i + params.st_nX) * params.dX + params.originX;
        for (index_t k = 0; k < params.phy_nZ; k++) {
            Real Z = real(k + params.st_nZ) * params.dZ + params.originZ;
            for (index_t j = 0; j < params.phy_nY; j++) {
                const Real Y = real(j + params.st_nY) * params.dY + params.originY;
                PU(interpData_ptr, front_i, j, k) = domData.westBF.UF(X, Y, Z, time);
                PV(interpData_ptr, front_i, j, k) = domData.westBF.VF(X, Y, Z, time);
                PW(interpData_ptr, front_i, j, k) = domData.westBF.WF(X, Y, Z, time);
                PP(interpData_ptr, front_i, j, k) = domData.westBF.PF(X, Y, Z, time);
            }
        }
    } else {
        for (index_t k = 0; k < params.phy_nZ; k++) {
            for (index_t j = 0; j < params.phy_nY; j++) {
                PU(interpData_ptr, front_i, j, k) = mu::interp_U_onGrid(data, front_i, j, k);
                PV(interpData_ptr, front_i, j, k) = mu::interp_V_onGrid(data, front_i, j, k);
                PW(interpData_ptr, front_i, j, k) = mu::interp_W_onGrid(data, front_i, j, k);
                PP(interpData_ptr, front_i, j, k) = mu::interp_P_onGrid(data, front_i, j, k);
            }
        }
    }

#define ITERATE_X_PHYSICAL_ROW(i, compute_physics)                                          \
    _Pragma("omp simd")                                                                     \
    ITERATE_ROW(i, 1, params.phy_nX - 1, compute_physics, X)

    // INTERPOLATE NORTH-EAST EDGE
    if (params.isOnTop || params.isOnRight) {
        auto &UF = params.isOnTop ? domData.northBF.UF : domData.eastBF.UF;
        auto &VF = params.isOnTop ? domData.northBF.VF : domData.eastBF.VF;
        auto &WF = params.isOnTop ? domData.northBF.WF : domData.eastBF.WF;
        auto &PF = params.isOnTop ? domData.northBF.PF : domData.eastBF.PF;

        const Real Y = real(north_j + params.st_nY) * params.dY + params.originY;
        const Real Z = real(east_k + params.st_nZ) * params.dZ + params.originZ;
        ITERATE_X_PHYSICAL_ROW(i, true)
            PU(interpData_ptr, i, north_j, east_k) = UF(X, Y, Z, time);
            PV(interpData_ptr, i, north_j, east_k) = VF(X, Y, Z, time);
            PW(interpData_ptr, i, north_j, east_k) = WF(X, Y, Z, time);
            PP(interpData_ptr, i, north_j, east_k) = PF(X, Y, Z, time);
        ITERATE_ROW_END()
    } else {
        ITERATE_X_PHYSICAL_ROW(i, false)
            PU(interpData_ptr, i, north_j, east_k) = mu::interp_U_onGrid(data, i, north_j, east_k);
            PV(interpData_ptr, i, north_j, east_k) = mu::interp_V_onGrid(data, i, north_j, east_k);
            PW(interpData_ptr, i, north_j, east_k) = mu::interp_W_onGrid(data, i, north_j, east_k);
            PP(interpData_ptr, i, north_j, east_k) = mu::interp_P_onGrid(data, i, north_j, east_k);
        ITERATE_ROW_END()
    }

    // INTERPOLATE NORTH-WEST EDGE
    if (params.isOnTop || params.isOnLeft) {
        auto &UF = params.isOnTop ? domData.northBF.UF : domData.westBF.UF;
        auto &VF = params.isOnTop ? domData.northBF.VF : domData.westBF.VF;
        auto &WF = params.isOnTop ? domData.northBF.WF : domData.westBF.WF;
        auto &PF = params.isOnTop ? domData.northBF.PF : domData.westBF.PF;

        const Real Y = real(north_j + params.st_nY) * params.dY + params.originY;
        const Real Z = real(west_k + params.st_nZ) * params.dZ + params.originZ;
        ITERATE_X_PHYSICAL_ROW(i, true)
            PU(interpData_ptr, i, north_j, west_k) = UF(X, Y, Z, time);
            PV(interpData_ptr, i, north_j, west_k) = VF(X, Y, Z, time);
            PW(interpData_ptr, i, north_j, west_k) = WF(X, Y, Z, time);
            PP(interpData_ptr, i, north_j, west_k) = PF(X, Y, Z, time);
        ITERATE_ROW_END()
    } else {
        ITERATE_X_PHYSICAL_ROW(i, false)
            PU(interpData_ptr, i, north_j, west_k) = mu::interp_U_onGrid(data, i, north_j, west_k);
            PV(interpData_ptr, i, north_j, west_k) = mu::interp_V_onGrid(data, i, north_j, west_k);
            PW(interpData_ptr, i, north_j, west_k) = mu::interp_W_onGrid(data, i, north_j, west_k);
            PP(interpData_ptr, i, north_j, west_k) = mu::interp_P_onGrid(data, i, north_j, west_k);
        ITERATE_ROW_END()
    }

    // INTERPOLATE SOUTH-EAST EDGE
    if (params.isOnBottom || params.isOnRight) {
        auto &UF = params.isOnBottom ? domData.southBF.UF : domData.eastBF.UF;
        auto &VF = params.isOnBottom ? domData.southBF.VF : domData.eastBF.VF;
        auto &WF = params.isOnBottom ? domData.southBF.WF : domData.eastBF.WF;
        auto &PF = params.isOnBottom ? domData.southBF.PF : domData.eastBF.PF;

        const Real Y = real(south_j + params.st_nY) * params.dY + params.originY;
        const Real Z = real(east_k + params.st_nZ) * params.dZ + params.originZ;
        ITERATE_X_PHYSICAL_ROW(i, true)
            PU(interpData_ptr, i, south_j, east_k) = UF(X, Y, Z, time);
            PV(interpData_ptr, i, south_j, east_k) = VF(X, Y, Z, time);
            PW(interpData_ptr, i, south_j, east_k) = WF(X, Y, Z, time);
            PP(interpData_ptr, i, south_j, east_k) = PF(X, Y, Z, time);
        ITERATE_ROW_END()
    } else {
        ITERATE_X_PHYSICAL_ROW(i, false)
            PU(interpData_ptr, i, south_j, east_k) = mu::interp_U_onGrid(data, i, south_j, east_k);
            PV(interpData_ptr, i, south_j, east_k) = mu::interp_V_onGrid(data, i, south_j, east_k);
            PW(interpData_ptr, i, south_j, east_k) = mu::interp_W_onGrid(data, i, south_j, east_k);
            PP(interpData_ptr, i, south_j, east_k) = mu::interp_P_onGrid(data, i, south_j, east_k);
        ITERATE_ROW_END()
    }

    // INTERPOLATE SOUTH-WEST EDGE
    if (params.isOnBottom || params.isOnLeft) {
        auto &UF = params.isOnBottom ? domData.southBF.UF : domData.westBF.UF;
        auto &VF = params.isOnBottom ? domData.southBF.VF : domData.westBF.VF;
        auto &WF = params.isOnBottom ? domData.southBF.WF : domData.westBF.WF;
        auto &PF = params.isOnBottom ? domData.southBF.PF : domData.westBF.PF;

        const Real Y = real(south_j + params.st_nY) * params.dY + params.originY;
        const Real Z = real(west_k + params.st_nZ) * params.dZ + params.originZ;
        ITERATE_X_PHYSICAL_ROW(i, true)
            PU(interpData_ptr, i, south_j, west_k) = UF(X, Y, Z, time);
            PV(interpData_ptr, i, south_j, west_k) = VF(X, Y, Z, time);
            PW(interpData_ptr, i, south_j, west_k) = WF(X, Y, Z, time);
            PP(interpData_ptr, i, south_j, west_k) = PF(X, Y, Z, time);
        ITERATE_ROW_END()
    } else {
        ITERATE_X_PHYSICAL_ROW(i, false)
            PU(interpData_ptr, i, south_j, west_k) = mu::interp_U_onGrid(data, i, south_j, west_k);
            PV(interpData_ptr, i, south_j, west_k) = mu::interp_V_onGrid(data, i, south_j, west_k);
            PW(interpData_ptr, i, south_j, west_k) = mu::interp_W_onGrid(data, i, south_j, west_k);
            PP(interpData_ptr, i, south_j, west_k) = mu::interp_P_onGrid(data, i, south_j, west_k);
        ITERATE_ROW_END()
    }
#undef ITERATE_X_PHYSICAL_ROW

}
#endif //INTERPOLATION_HPP
