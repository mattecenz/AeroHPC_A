#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include "SolverData.hpp"
#include "ForcingTerm.hpp"
#include "Boundaries.hpp"
#include "mathUtils.hpp"
#include "PressureSolver.hpp"

#include "printBuffer.hpp"

namespace mu = mathUtils;

////RHS function
#define compute_rhs(C)                                                                                            \
inline Real compute_rhs_##C(Real *data, const Real nu, const index_t i, const index_t j, const index_t k) \
{                                                                                                     \
return -mu::conv_##C(data, i, j, k) + nu * mu::lap_##C(data, i, j, k);                            \
}

compute_rhs(U)

compute_rhs(V)

compute_rhs(W)


/// FORCING TERM ///////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef ForcingT
#define getPhys(i, j, k)                         \
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
#define ITERATE_DOMAIN_VELOCITY(i,j,k)                                      \
const index_t nz = params.isOnRight ? params.loc_nZ - 1 : params.loc_nZ;    \
const index_t ny = params.isOnTop ? params.loc_nY - 1 : params.loc_nY;      \
const index_t nx = params.loc_nY - 1;                                       \
for (index_t k = 0; k < nz; ++k) {                                          \
    for (index_t j = 0; j < ny; ++j) {                                      \
    _Pragma("omp simd")                                                     \
        for (index_t i = 0; i < nx; ++i) {

#define ITERATE_DOMAIN_PRESSURE(i,j,k)                                  \
const index_t nz = params.loc_nZ;                                       \
const index_t ny = params.loc_nY;                                       \
const index_t nx = params.loc_nX;                                       \
for (index_t k = 0; k < nz; ++k) {                                      \
    for (index_t j = 0; j < ny; ++j) {                                  \
    _Pragma("omp simd")                                                 \
        for (index_t i = 0; i < nx; ++i) {

#define ITERATE_DOMAIN_END()                                        \
        }                                                           \
    }                                                               \
}


#define Y2star_C(C, Y2star, U_N, P_N)                   \
ITERATE_DOMAIN_VELOCITY(i, j, k)                        \
    getForce##C(force, i, j, k);                        \
    getPressureGrad##C(d_press, P_N, i, j, k);          \
    const Real r = compute_rhs_##C(U_N, nu, i, j, k);   \
    rhs_##C(i, j, k) = (r + force);                     \
    C(Y2star, i, j, k) = C(U_N, i, j, k)                \
                            + k_0 * (r + force)         \
                            - k_0 * d_press;            \
ITERATE_DOMAIN_END()

#define Y2star(Y2star, U_N, P_N)    \
    {Y2star_C(U, Y2star, U_N, P_N)} \
    {Y2star_C(V, Y2star, U_N, P_N)} \
    {Y2star_C(W, Y2star, U_N, P_N)}

#define Y2(Y2, Y2star, PHI_2, PHI_2_P_N, P_N)                                           \
{ITERATE_DOMAIN_VELOCITY(i, j, k)                                                       \
    U(Y2, i, j, k) = U(Y2star, i, j, k) - params.dt * mu::dp_dx_U(PHI_2_P_N, i, j, k);  \
    V(Y2, i, j, k) = V(Y2star, i, j, k) - params.dt * mu::dp_dy_V(PHI_2_P_N, i, j, k);  \
    W(Y2, i, j, k) = W(Y2star, i, j, k) - params.dt * mu::dp_dz_W(PHI_2_P_N, i, j, k);  \
ITERATE_DOMAIN_END()}                                                                   \
{ITERATE_DOMAIN_PRESSURE(i, j, k)                                                       \
    P(PHI_2, i, j, k) = P(P_N, i, j, k) + P(PHI_2_P_N, i, j, k);                        \
ITERATE_DOMAIN_END()}

#define Y3star_C(C, Y3star, Y2, PHI_2)                  \
ITERATE_DOMAIN_VELOCITY(i, j, k)                        \
    getForce##C(force, i, j, k);                        \
    getPressureGrad##C(d_press, PHI_2, i, j, k);        \
    const Real r1 = rhs_##C(i, j, k);                   \
    const Real r2 = compute_rhs_##C(Y2, nu, i, j, k);   \
    rhs_##C(i, j, k) = (r2 + force);                    \
    C(Y3star, i, j, k) = C(Y2, i, j, k)                 \
                       - k_1 * r1                       \
                       + k_2 * (r2 + force)             \
                       - k_3 * d_press;                 \
ITERATE_DOMAIN_END()

#define Y3star(Y3star, Y2, PHI_2)  \
    {Y3star_C(U, Y3star, Y2, PHI_2)} \
    {Y3star_C(V, Y3star, Y2, PHI_2)} \
    {Y3star_C(W, Y3star, Y2, PHI_2)}

#define Y3(Y3, Y3star, PHI_3, PHI_3_PHI_2, PHI_2)                                           \
{ITERATE_DOMAIN_VELOCITY(i, j, k)                                                           \
    U(Y3, i, j, k) = U(Y3star, i, j, k) - params.dt * mu::dp_dx_U(PHI_3_PHI_2, i, j, k);    \
    V(Y3, i, j, k) = V(Y3star, i, j, k) - params.dt * mu::dp_dy_V(PHI_3_PHI_2, i, j, k);    \
    W(Y3, i, j, k) = W(Y3star, i, j, k) - params.dt * mu::dp_dz_W(PHI_3_PHI_2, i, j, k);    \
ITERATE_DOMAIN_END()}                                                                       \
{ITERATE_DOMAIN_PRESSURE(i, j, k)                                                           \
    P(PHI_3, i, j, k) = P(PHI_2, i, j, k) + P(PHI_3_PHI_2, i, j, k);                        \
ITERATE_DOMAIN_END()}


#define U_N1star_C(C, U_N1star, Y3, PHI_3)              \
ITERATE_DOMAIN_VELOCITY(i, j, k)                        \
    getForce##C(force, i, j, k);                        \
    getPressureGrad##C(d_press, PHI_3, i, j, k);        \
    const Real r = compute_rhs_##C(Y3, nu, i, j, k);    \
    C(U_N1star, i, j, k) = C(Y3, i, j, k)               \
                            - k_4 * rhs_##C(i, j, k)    \
                            + k_5 * (r + force)         \
                            - k_6 * d_press;            \
ITERATE_DOMAIN_END()

#define U_N1star(U_N1star, Y3, PHI_3)       \
    {U_N1star_C(U, U_N1star, Y3, PHI_3)}    \
    {U_N1star_C(V, U_N1star, Y3, PHI_3)}    \
    {U_N1star_C(W, U_N1star, Y3, PHI_3)}

#define U_N1(U_N1, U_N1star, P_N1, P_N1_PHI_3, PHI_3)                                           \
{ITERATE_DOMAIN_VELOCITY(i, j, k)                                                               \
    U(U_N1, i, j, k) = U(U_N1star, i, j, k) - params.dt * mu::dp_dx_U(P_N1_PHI_3, i, j, k);     \
    V(U_N1, i, j, k) = V(U_N1star, i, j, k) - params.dt * mu::dp_dy_V(P_N1_PHI_3, i, j, k);     \
    W(U_N1, i, j, k) = W(U_N1star, i, j, k) - params.dt * mu::dp_dz_W(P_N1_PHI_3, i, j, k);     \
ITERATE_DOMAIN_END()}                                                                           \
{ITERATE_DOMAIN_PRESSURE(i, j, k)                                                               \
    P(P_N1, i, j, k) = P(PHI_3, i, j, k) + P(P_N1_PHI_3, i, j, k);                              \
ITERATE_DOMAIN_END()}


#define Load_B(constant, VELOCITY)                                  \
{ITERATE_DOMAIN_PRESSURE(i, j, k)                                   \
    rhs_P(i, j, k) = constant * mu::vel_div(VELOCITY, i, j, k);     \
ITERATE_DOMAIN_END()}

#define Unload_B(PRESSURE)                  \
{ITERATE_DOMAIN_PRESSURE(i, j, k)           \
    P(PRESSURE, i, j, k) = rhs_P(i, j, k);  \
ITERATE_DOMAIN_END()}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct RKConst {
    static constexpr Real alpha0 = 64.0 / 120.0;
    static constexpr Real alpha1 = 34.0 / 120.0;
    static constexpr Real alpha2 = 50.0 / 120.0;
    static constexpr Real alpha3 = alpha2 - alpha1;
    static constexpr Real alpha4 = 50.0 / 120.0;
    static constexpr Real alpha5 = 90.0 / 120.0;
    static constexpr Real alpha6 = alpha5 - alpha4;
    static constexpr Real beta0 = 64.0 / 120.0;
    static constexpr Real beta1 = 80.0 / 120.0;
};

// Runge-Kutta method
inline void rungeKutta(const Real time) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //b_print(model, dir, 0);

    const Real nu = (real(1) / params.Re);

    // kappa -> weighted_deltat
    const Real k_0 = RKConst::alpha0 * params.dt;
    const Real k_1 = RKConst::alpha1 * params.dt;
    const Real k_2 = RKConst::alpha2 * params.dt;
    const Real k_3 = RKConst::alpha3 * params.dt;
    const Real k_4 = RKConst::alpha4 * params.dt;
    const Real k_5 = RKConst::alpha5 * params.dt;
    const Real k_6 = RKConst::alpha6 * params.dt;
    const Real inv_k_0 = real(1.0) / k_0;
    const Real inv_k_3 = real(1.0) / k_3;
    const Real inv_k_6 = real(1.0) / k_6;
    const Real t_0 = time + RKConst::beta0 * params.dt;
    const Real t_1 = time + RKConst::beta1 * params.dt;
    const Real t_2 = time + params.dt;


#ifdef ForcingT
    ForcingTerm ft(params.Re, time);
#endif

    /// Y2* //////////////////////////////////////////////////////////////////////////////////////////////
    {
        Y2star(rkData.buffer_data, rkData.model_data, rkData.model_data)

        //b_print(model_buff, 1);
        // TODO APPLY BOUNDARIES ON VELOCITY
        //boundary_cond.apply(model_buff, t_0);
        //b_print(model_buff, dir, 2);
    }

#ifndef DISABLE_PRESSURE
    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        // TODO LOAD PRESSURE BUFFER
        Load_B(params.dt, rkData.buffer_data)

        // SOLVE FOR phi2-pn
        // TODO CALL PRESSURE SOLVER
        solvePressure();

        // TODO UNLOAD PRESSURE BUFFER
        Unload_B(rkData.buffer_data)

        //b_print(model_buff, dir, 3);
        // TODO APPLY BOUNDARY ON PRESSURE
        // boundary_cond.apply(model_buff, t_0);
        //b_print(model_buff, dir, 4);
    }


    /// Y2 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        Y2(rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.model_data)
        // TODO ADD EXPLICIT CALL FOR PHI2

        //b_print(model_buff, dir, 5);
        // TODO APPLY BOUNDARY ON VELOCITY
        // TODO APPLY BOUNDARY ON PRESSURE
        // boundary_cond.apply(model_buff, t_0);
        //b_print(model_buff, dir, 6);
    }
#endif

#ifdef ForcingT
    ft.set_time(t_0);
#endif

    /// Y3* //////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        Y3star(rkData.model_data, rkData.buffer_data, rkData.buffer_data)

        //TODO boundary of velocity
        //b_print(model, dir, 7);
        //boundary_cond.apply(model, t_1);
        //b_print(model, dir, 8);
    }

#ifndef DISABLE_PRESSURE
    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        Load_B(params.dt, rkData.model_data)

        solvePressure();

        Unload_B(rhs_buff, model)

        //TODO boundary of pressure
        //b_print(model, dir, 9);
        //boundary_cond.apply(model, t_1);
        //b_print(model, dir, 10);
    }

    /// Y3 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        Y3(rkData.model_data, rkData.model_data, rkData.model_data, rkData.model_data, rkData.buffer_data)

        //TODO boundary of velocity
        //TODO boundary of pressure
        //b_print(model, dir, 11);
        //boundary_cond.apply(model, t_1);
        //b_print(model, dir, 12);
    }
#endif


#ifdef ForcingT
    ft.set_time(t_1);
#endif

    /// u(n+1)* //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        U_N1star(rkData.buffer_data, rkData.model_data, rkData.model_data)

        //TODO boundary of velocity
        //b_print(model_buff, dir, 13);
        //boundary_cond.apply(model_buff, t_2);
        //b_print(model_buff, dir, 14);
    }

#ifndef DISABLE_PRESSURE
    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        Load_B(params.dt, rkData.buffer_data)

        solvePressure();

        Unload_B(rkData.buffer_data)

        //TODO boundary of pressure
        //b_print(model_buff, dir, 15);
        //boundary_cond.apply(model_buff, t_2);
        //b_print(model_buff, dir, 16);
    }

    /// un+1 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        U_N1(rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.model_data);

        //TODO boundary of velocity
        //TODO boundary of pressure
        //b_print(model_buff, dir, 17);
        //boundary_cond.apply(model_buff, t_2);
        //b_print(model_buff, dir, 18);
    }
#endif

    swap(rkData.buffer_data, rkData.model_data);
}

#endif // AEROHPC_A_RUNGEKUTTA_H
