#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include "GridData.hpp"
#include "ForcingTerm.hpp"
#include "Boundaries.hpp"
#include "mathUtils.hpp"
#include "poissonSolver.hpp"


/// TODO TO REMOVE ONLY FOR DEBUG ///////////////////////////////////////////////////
#include "printBuffer.hpp"
#ifdef DEBUG_PRINT_BUFFERS
#define b_print(buff, dir, n) \
if (!rank) { \
    std::string nn{dir + to_string(n)}; \
    print(buff, nn); \
}
#define c_dir(dir) \
if (!rank) {\
    create_directories(dir); \
}
#else
#define b_print(a,b,c) //
#define c_dir(a) //
#endif
/////////////////////////////////////////////////////////////////////////////////////


namespace mu = mathUtils;

////RHS function
#define rhs(C)                                                                                            \
inline Real rhs_##C(GridData &grid, const Real nu, const index_t i, const index_t j, const index_t k) \
{                                                                                                     \
return -mu::conv_##C(grid, i, j, k) + nu * mu::lap_##C(grid, i, j, k);                            \
}

rhs(U)

rhs(V)

rhs(W)


/// FORCING TERM ///////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef ForcingT
#define getPhys(i, j, k)                         \
    Real px = real(i + model.structure.px) * dx; \
    Real py = real(j + model.structure.py) * dy; \
    Real pz = real(k + model.structure.pz) * dz
#define getForceU(force, i, j, k) getPhys(i, j, k); const Real force = ft.computeGx(px + dx, py + sdy, pz + sdz)
#define getForceV(force, i, j, k) getPhys(i, j, k); const Real force = ft.computeGy(px + sdx, py + dy, pz + sdz)
#define getForceW(force, i, j, k) getPhys(i, j, k); const Real force = ft.computeGz(px + sdx, py + sdy, pz + dz)
#else
#define getForceU(force, i, j, k) constexpr Real force = 0
#define getForceV(force, i, j, k) constexpr Real force = 0
#define getForceW(force, i, j, k) constexpr Real force = 0
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// PRESSURE TERM //////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_PRESSURE
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
#define ITERATE_OVER_ALL_POINTS_START(i,j,k)    \
for (index_t k = 0; k < nz; ++k) {              \
    for (index_t j = 0; j < ny; ++j) {          \
        for (index_t i = 0; i < nx; ++i) {
#define ITERATE_OVER_ALL_POINTS_END()           \
        }                                       \
    }                                           \
}

#define Y2star(C, Y2star, U_N)                      \
ITERATE_OVER_ALL_POINTS_START(i, j, k)              \
    getForce##C(force, i, j, k);                    \
    getPressureGrad##C(d_press, U_N, i, j, k);      \
    const Real r = rhs_##C(U_N, nu, i, j, k);       \
    rhs_buff.C(i, j, k) = (r + force);              \
    Y2star.C(i, j, k) = U_N.C(i, j, k)              \
                            + k_0 * (r + force)     \
                            - k_0 * d_press;        \
ITERATE_OVER_ALL_POINTS_END()

#define Y2(Y2, Y2star, U_N)                                                 \
ITERATE_OVER_ALL_POINTS_START(i, j, k)                                      \
    Y2.U(i, j, k) = Y2star.U(i, j, k) - k_0 * mu::dp_dx_U(Y2star, i, j, k); \
    Y2.V(i, j, k) = Y2star.V(i, j, k) - k_0 * mu::dp_dy_V(Y2star, i, j, k); \
    Y2.W(i, j, k) = Y2star.W(i, j, k) - k_0 * mu::dp_dz_W(Y2star, i, j, k); \
ITERATE_OVER_ALL_POINTS_END()                                               \
ITERATE_OVER_ALL_POINTS_START(i, j, k)                                      \
    Y2.P(i, j, k) = Y2star.P(i, j, k) + U_N.P(i, j, k);                     \
ITERATE_OVER_ALL_POINTS_END()

#define Y3star(C, Y3star, Y2)                           \
ITERATE_OVER_ALL_POINTS_START(i, j, k)                  \
    getForce##C(force, i, j, k);                        \
    getPressureGrad##C(d_press, Y2, i, j, k);           \
    const Real r1 = rhs_buff.C(i, j, k);                \
    const Real r2 = rhs_##C(Y2, nu, i, j, k);           \
    rhs_buff.C(i, j, k) = (r2 + force);                 \
    Y3star.C(i, j, k) = Y2.C(i, j, k)                   \
                       - k_1 * r1                       \
                       + k_2 * (r2 + force)             \
                       - k_3 * d_press;                 \
ITERATE_OVER_ALL_POINTS_END()

#define Y3(Y3, Y3star, Y2)                                                  \
ITERATE_OVER_ALL_POINTS_START(i, j, k)                                      \
    Y3.U(i, j, k) = Y3star.U(i, j, k) - k_3 * mu::dp_dx_U(Y3star, i, j, k); \
    Y3.V(i, j, k) = Y3star.V(i, j, k) - k_3 * mu::dp_dy_V(Y3star, i, j, k); \
    Y3.W(i, j, k) = Y3star.W(i, j, k) - k_3 * mu::dp_dz_W(Y3star, i, j, k); \
ITERATE_OVER_ALL_POINTS_END()                                               \
ITERATE_OVER_ALL_POINTS_START(i, j, k)                                      \
    Y3.P(i, j, k) = Y3star.P(i, j, k) + Y2.P(i, j, k);                      \
ITERATE_OVER_ALL_POINTS_END()


#define U_N1star(C, U_N1, Y3)                           \
ITERATE_OVER_ALL_POINTS_START(i, j, k)                  \
    getForce##C(force, i, j, k);                        \
    getPressureGrad##C(d_press, Y3, i, j, k);           \
    const Real r = rhs_##C(Y3, nu, i, j, k);            \
    U_N1.C(i, j, k) = Y3.C(i, j, k)                     \
                            - k_4 * rhs_buff.C(i, j, k) \
                            + k_5 * (r + force)         \
                            - k_6 * d_press;            \
ITERATE_OVER_ALL_POINTS_END()

#define U_N1(U_N1, U_N1star, Y3)                                                    \
ITERATE_OVER_ALL_POINTS_START(i, j, k)                                              \
    U_N1.U(i, j, k) = U_N1star.U(i, j, k) - k_6 * mu::dp_dx_U(U_N1star, i, j, k);   \
    U_N1.V(i, j, k) = U_N1star.V(i, j, k) - k_6 * mu::dp_dy_V(U_N1star, i, j, k);   \
    U_N1.W(i, j, k) = U_N1star.W(i, j, k) - k_6 * mu::dp_dz_W(U_N1star, i, j, k);   \
ITERATE_OVER_ALL_POINTS_END()                                                       \
ITERATE_OVER_ALL_POINTS_START(i, j, k)                                              \
    U_N1.P(i, j, k) = U_N1star.P(i, j, k) + Y3.P(i, j, k);                          \
ITERATE_OVER_ALL_POINTS_END()


#define Load_B(constant, b_buffer, velocity_in)                                     \
ITERATE_OVER_ALL_POINTS_START(i, j, k)                                              \
    b_buffer.P(i, j, k) = constant * mu::vel_div(velocity_in, i, j, k);             \
ITERATE_OVER_ALL_POINTS_END()

#define Unload_B(b_buffer, pressure_out)            \
ITERATE_OVER_ALL_POINTS_START(i, j, k)              \
    pressure_out.P(i, j, k) = b_buffer.P(i, j, k);  \
ITERATE_OVER_ALL_POINTS_END()
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
inline void rungeKutta(GridData &model, GridData &model_buff, GridData &rhs_buff,
                       Real reynolds, Real deltat, index_t iteration,
                       Boundaries &boundary_cond, poissonSolver &p_solver) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string dir = "./iterations/" + to_string(iteration) + "/";

    c_dir(dir);
    b_print(model, dir, 0);

    const Real time = deltat * real(iteration);

    const Real nu = (real(1) / reynolds);

    const Real dx = model.structure.dx;
    const Real dy = model.structure.dy;
    const Real dz = model.structure.dz;

    const Real sdx = model.structure.sdx;
    const Real sdy = model.structure.sdy;
    const Real sdz = model.structure.sdz;

    const index_t nx = model.structure.nx;
    const index_t ny = model.structure.ny;
    const index_t nz = model.structure.nz;

    // kappa -> weighted_deltat
    const Real k_0 = RKConst::alpha0 * deltat;
    const Real k_1 = RKConst::alpha1 * deltat;
    const Real k_2 = RKConst::alpha2 * deltat;
    const Real k_3 = RKConst::alpha3 * deltat;
    const Real k_4 = RKConst::alpha4 * deltat;
    const Real k_5 = RKConst::alpha5 * deltat;
    const Real k_6 = RKConst::alpha6 * deltat;
    const Real inv_k_0 = real(1.0) / k_0;
    const Real inv_k_3 = real(1.0) / k_3;
    const Real inv_k_6 = real(1.0) / k_6;
    const Real t_0 = time + RKConst::beta0 * deltat;
    const Real t_1 = time + RKConst::beta1 * deltat;
    const Real t_2 = time + deltat;


#ifdef ForcingT
    ForcingTerm ft(reynolds, time);
#endif

    /// Y2* //////////////////////////////////////////////////////////////////////////////////////////////
    {
        Y2star(U, model_buff, model);
        Y2star(V, model_buff, model);
        Y2star(W, model_buff, model);

        b_print(model_buff, dir, 1);
        boundary_cond.apply(model_buff, t_0);
        b_print(model_buff, dir, 2);
    }

#ifdef ENABLE_PRESSURE
    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        Load_B(inv_k_0, rhs_buff, model_buff)

        // SOLVE FOR phi2-pn
        p_solver.solve(&rhs_buff.P(0,0,0));
    }

    /// UPDATE PRESSURE /////////////////////////////////////////////////////////////////////////////////////////////
    {
        Unload_B(rhs_buff, model_buff)

        b_print(model_buff, dir, 3);
        boundary_cond.apply(model_buff, t_0);
        b_print(model_buff, dir, 4);
    }


    /// Y2 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        Y2(model_buff, model_buff, model);

        b_print(model_buff, dir, 5);
        boundary_cond.apply(model_buff, t_0);
        b_print(model_buff, dir, 6);
    }
#endif

#ifdef ForcingT
    ft.set_time(t_0);
#endif

    /// Y3* //////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        Y3star(U, model, model_buff);
        Y3star(V, model, model_buff);
        Y3star(W, model, model_buff);

        b_print(model, dir, 7);
        boundary_cond.apply(model, t_1);
        b_print(model, dir, 8);
    }

#ifdef ENABLE_PRESSURE
    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        Load_B(inv_k_3, rhs_buff, model)

        p_solver.solve(&rhs_buff.P(0,0,0));
    }

    /// UPDATE PRESSURE /////////////////////////////////////////////////////////////////////////////////////////////
    {
        Unload_B(rhs_buff, model)

        b_print(model, dir, 9);
        boundary_cond.apply(model, t_1);
        b_print(model, dir, 10);
    }

    /// Y3 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        Y3(model, model, model_buff);

        b_print(model, dir, 11);
        boundary_cond.apply(model, t_1);
        b_print(model, dir, 12);
    }
#endif


#ifdef ForcingT
    ft.set_time(t_1);
#endif

    /// u(n+1)* //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        U_N1star(U, model_buff, model);
        U_N1star(V, model_buff, model);
        U_N1star(W, model_buff, model);


        b_print(model_buff, dir, 13);
        boundary_cond.apply(model_buff, t_2);
        b_print(model_buff, dir, 14);
    }

#ifdef ENABLE_PRESSURE
    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        Load_B(inv_k_6, rhs_buff, model_buff)

        p_solver.solve(&rhs_buff.P(0,0,0));
    }

    /// UPDATE PRESSURE /////////////////////////////////////////////////////////////////////////////////////////////
    {
        Unload_B(rhs_buff, model_buff)

        b_print(model_buff, dir, 15);
        boundary_cond.apply(model_buff, t_2);
        b_print(model_buff, dir, 16);
    }

    /// un+1 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        U_N1(model_buff, model_buff, model);

        b_print(model_buff, dir, 17);
        boundary_cond.apply(model_buff, t_2);
        b_print(model_buff, dir, 18);
    }
#endif

    model.swap(model_buff);
}

#endif // AEROHPC_A_RUNGEKUTTA_H
