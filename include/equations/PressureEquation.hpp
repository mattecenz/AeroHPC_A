#ifndef PRESSUREEQUATION_HPP
#define PRESSUREEQUATION_HPP

#include "PressureSolver.hpp"
#include "utils/macroUtils.hpp"

/// PRESSURE EQUATIONS /////////////////////////////////////////////////////////////////////////////
#define Y2(Y2, Y2star, PHI_2, PHI_2_P_N, P_N)                                           \
ITERATE_DOMAIN_VELOCITY(i, j, k, false, NO_SKIP)                                        \
    U(Y2, i, j, k) = U(Y2star, i, j, k) - 1.0 * mu::dp_dx_U(PHI_2_P_N, i, j, k);  \
    V(Y2, i, j, k) = V(Y2star, i, j, k) - 1.0 * mu::dp_dy_V(PHI_2_P_N, i, j, k);  \
    W(Y2, i, j, k) = W(Y2star, i, j, k) - 1.0 * mu::dp_dz_W(PHI_2_P_N, i, j, k);  \
ITERATE_DOMAIN_END()                                                                    \
ITERATE_DOMAIN_PRESSURE(i, j, k, false)                                                 \
    P(PHI_2, i, j, k) = P(P_N, i, j, k) + P(PHI_2_P_N, i, j, k);                        \
ITERATE_DOMAIN_END()


#define Y3(Y3, Y3star, PHI_3, PHI_3_PHI_2, PHI_2)                                           \
ITERATE_DOMAIN_VELOCITY(i, j, k, false, NO_SKIP)                                            \
    U(Y3, i, j, k) = U(Y3star, i, j, k) - 1.0 * mu::dp_dx_U(PHI_3_PHI_2, i, j, k);    \
    V(Y3, i, j, k) = V(Y3star, i, j, k) - 1.0 * mu::dp_dy_V(PHI_3_PHI_2, i, j, k);    \
    W(Y3, i, j, k) = W(Y3star, i, j, k) - 1.0 * mu::dp_dz_W(PHI_3_PHI_2, i, j, k);    \
ITERATE_DOMAIN_END()                                                                        \
ITERATE_DOMAIN_PRESSURE(i, j, k, false)                                                     \
    P(PHI_3, i, j, k) = P(PHI_2, i, j, k) + P(PHI_3_PHI_2, i, j, k);                        \
ITERATE_DOMAIN_END()


#define U_N1(U_N1, U_N1star, P_N1, P_N1_PHI_3, PHI_3)                                           \
ITERATE_DOMAIN_VELOCITY(i, j, k, false, NO_SKIP)                                                \
    U(U_N1, i, j, k) = U(U_N1star, i, j, k) - 1.0 * mu::dp_dx_U(P_N1_PHI_3, i, j, k);     \
    V(U_N1, i, j, k) = V(U_N1star, i, j, k) - 1.0 * mu::dp_dy_V(P_N1_PHI_3, i, j, k);     \
    W(U_N1, i, j, k) = W(U_N1star, i, j, k) - 1.0 * mu::dp_dz_W(P_N1_PHI_3, i, j, k);     \
ITERATE_DOMAIN_END()                                                                            \
ITERATE_DOMAIN_PRESSURE(i, j, k, false)                                                         \
    P(P_N1, i, j, k) = P(PHI_3, i, j, k) + P(P_N1_PHI_3, i, j, k);                              \
ITERATE_DOMAIN_END()


#define Load_B(constant, VELOCITY)                                  \
ITERATE_DOMAIN_PRESSURE(i, j, k, false)                             \
    rhs_P(k, j, i) = constant * mu::vel_div(VELOCITY, i, j, k);     \
ITERATE_DOMAIN_END()


#define Unload_B(PRESSURE)                  \
ITERATE_DOMAIN_PRESSURE(i, j, k, false)     \
    P(PRESSURE, k, j, i) = rhs_P(i, j, k);  \
ITERATE_DOMAIN_END()
//////////////////////////////////////////////////////////////////////////////////////////////////////


/// PRESSURE SOLVER MACRO ////////////////////////////////////////////////////////////////////////////
#define P_Eq(constant, VELOCITY_IN, PRESSURE_OUT)   \
    Load_B(constant, VELOCITY_IN)                   \
    solvePressure();                                \
    Unload_B(PRESSURE_OUT)
//////////////////////////////////////////////////////////////////////////////////////////////////////

#endif //PRESSUREEQUATION_HPP
