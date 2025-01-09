#ifndef MOMENTUMEQUATION_HPP
#define MOMENTUMEQUATION_HPP

#include "Traits.hpp"
#include "data/SolverData.hpp"
#include "utils/mathUtils.hpp"
#include "utils/macroUtils.hpp"

namespace mu = mathUtils;

////RHS function
#define compute_rhs(C)                                                                                          \
    inline Real compute_rhs_##C(Real *data, const Real nu, const index_t i, const index_t j, const index_t k)   \
    {                                                                                                           \
        return -mu::conv_##C(data, i, j, k) + consts.nu * mu::lap_##C(data, i, j, k);                           \
    }

compute_rhs(U)

compute_rhs(V)

compute_rhs(W)


/// FORCING TERM ///////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef ForcingT
#define getPhys(i, j, k)                          \
    Real px = real(i + params.st_nX) * params.dX; \
    Real py = real(j + params.st_nY) * params.dY; \
    Real pz = real(k + params.st_nZ) * params.dZ
#define getForceU(force, i, j, k) getPhys(i, j, k); const Real force = ft.computeGx(px + params.dX, py + params.dX2, pz + params.dZ2)
#define getForceV(force, i, j, k) getPhys(i, j, k); const Real force = ft.computeGy(px + params.dX2, py + params.dY, pz + params.dZ2)
#define getForceW(force, i, j, k) getPhys(i, j, k); const Real force = ft.computeGz(px + params.dX2, py + params.dY2, pz + params.dZ)
#else
#define getForceU(force, i, j, k) constexpr Real force = 0
#define getForceV(force, i, j, k) constexpr Real force = 0
#define getForceW(force, i, j, k) constexpr Real force = 0
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// PRESSURE TERM //////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef DISABLE_PRESSURE
#define getPressureGradU(d_press, buff, i, j, k) const Real d_press = mu::dp_dx_U(buff, i, j, k)
#define getPressureGradV(d_press, buff, i, j, k) const Real d_press = mu::dp_dy_V(buff, i, j, k)
#define getPressureGradW(d_press, buff, i, j, k) const Real d_press = mu::dp_dz_W(buff, i, j, k)
#else
#define getPressureGradU(d_press, buff, i, j, k) constexpr Real d_press = 0
#define getPressureGradV(d_press, buff, i, j, k) constexpr Real d_press = 0
#define getPressureGradW(d_press, buff, i, j, k) constexpr Real d_press = 0
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// RK STEPS ///////////////////////////////////////////////////////////////////////////////////////////////////////////
#define Y2star_C(C, Y2star, U_N, P_N)                           \
ITERATE_DOMAIN_VELOCITY(i, j, k)                                \
    getForce##C(force, i, j, k);                                \
    getPressureGrad##C(d_press, P_N, i, j, k);                  \
    const Real r = compute_rhs_##C(U_N, consts.nu, i, j, k);    \
    rhs_##C(i, j, k) = (r + force);                             \
    C(Y2star, i, j, k) = C(U_N, i, j, k)                        \
                            + consts.k_0 * (r + force)          \
                            - consts.k_0 * d_press;             \
ITERATE_DOMAIN_END()

#define Y2star(Y2star, U_N, P_N)    \
    {Y2star_C(U, Y2star, U_N, P_N)} \
    {Y2star_C(V, Y2star, U_N, P_N)} \
    {Y2star_C(W, Y2star, U_N, P_N)}



#define Y3star_C(C, Y3star, Y2, PHI_2)                  \
ITERATE_DOMAIN_VELOCITY(i, j, k)                        \
    getForce##C(force, i, j, k);                        \
    getPressureGrad##C(d_press, PHI_2, i, j, k);        \
    const Real r1 = rhs_##C(i, j, k);                   \
    const Real r2 = compute_rhs_##C(Y2, consts.nu, i, j, k);   \
    rhs_##C(i, j, k) = (r2 + force);                    \
    C(Y3star, i, j, k) = C(Y2, i, j, k)                 \
                       - consts.k_1 * r1                       \
                       + consts.k_2 * (r2 + force)             \
                       - consts.k_3 * d_press;                 \
ITERATE_DOMAIN_END()

#define Y3star(Y3star, Y2, PHI_2)  \
    {Y3star_C(U, Y3star, Y2, PHI_2)} \
    {Y3star_C(V, Y3star, Y2, PHI_2)} \
    {Y3star_C(W, Y3star, Y2, PHI_2)}




#define U_N1star_C(C, U_N1star, Y3, PHI_3)              \
ITERATE_DOMAIN_VELOCITY(i, j, k)                        \
    getForce##C(force, i, j, k);                        \
    getPressureGrad##C(d_press, PHI_3, i, j, k);        \
    const Real r = compute_rhs_##C(Y3, consts.nu, i, j, k);    \
    C(U_N1star, i, j, k) = C(Y3, i, j, k)               \
                            - consts.k_4 * rhs_##C(i, j, k)    \
                            + consts.k_5 * (r + force)         \
                            - consts.k_6 * d_press;            \
ITERATE_DOMAIN_END()

#define U_N1star(U_N1star, Y3, PHI_3)       \
    {U_N1star_C(U, U_N1star, Y3, PHI_3)}    \
    {U_N1star_C(V, U_N1star, Y3, PHI_3)}    \
    {U_N1star_C(W, U_N1star, Y3, PHI_3)}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif //MOMENTUMEQUATION_HPP
