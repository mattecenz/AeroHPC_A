#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include "GridData.hpp"
#include "ForcingTerm.hpp"
#include "Boundaries.hpp"
#include "mathUtils.hpp"

#ifdef ForcingT
#define getPhys(i, j, k)                         \
    Real px = real(i + model.structure.px) * dx; \
    Real py = real(j + model.structure.py) * dy; \
    Real pz = real(k + model.structure.pz) * dz

#define getForceU(ftt) ftt.computeGx(px + dx, py + sdy, pz + sdz)

#define getForceV(ftt) ftt.computeGy(px + sdx, py + dy, pz + sdz)

#define getForceW(ftt) ftt.computeGz(px + sdx, py + sdy, pz + dz)
#endif

namespace mu = mathUtils;

struct RKConst
{
    static constexpr Real alpha0 = 8.0 / 15.0;
    static constexpr Real alpha1 = 17.0 / 60.0;
    static constexpr Real alpha2 = 5.0 / 12.0;
    static constexpr Real alpha3 = 3.0 / 4.0;
    static constexpr Real alpha4 = 2.0 / 3.0;
    static constexpr Real alpha5 = 2.0 / 15.0;
    static constexpr Real alpha6 = 1.0 / 3.0;
    static constexpr Real alpha7 = 120.0 / 64.0;
    static constexpr Real alpha8 = 120.0 / 16.0;
    static constexpr Real alpha9 = 120.0 / 40.0;
};

////RHS function
#define rhs(C)                                                                                            \
    inline Real rhs_##C(GridData &grid, const Real nu, const index_t i, const index_t j, const index_t k) \
    {                                                                                                     \
        return -mu::conv_##C(grid, i, j, k) + nu * mu::lap_##C(grid, i, j, k);                            \
    }

rhs(U)

    rhs(V)

        rhs(W)

    // Runge-Kutta method
    void rungeKutta(GridData &model, GridData &model_buff, GridData &rhs_buff,
                    Real reynolds, Real deltat, Real time,
                    Boundaries &boundary_cond)
{

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
    std::array<const Real, 10> kappa{
        RKConst::alpha0 * deltat,
        RKConst::alpha1 * deltat,
        RKConst::alpha2 * deltat,
        RKConst::alpha3 * deltat,
        RKConst::alpha4 * deltat,
        RKConst::alpha5 * deltat,
        RKConst::alpha6 * deltat,
        RKConst::alpha7 * deltat,
        RKConst::alpha8 * deltat,
        RKConst::alpha9 * deltat};

#ifdef ForcingT
    ForcingTerm ft(reynolds, time);
#endif
    // What parameters (the first 2 are reals related to the grid size, but ours is different)
    // poissonSolver solver();
    double X[nx * ny * nz];
    double b[nx * ny * nz];

    /// Y2* //////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceU(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_U(model, nu, i, j, k);
                    rhs_buff.U(i, j, k) = (r + force);

                    model_buff.U(i, j, k) = model.U(i, j, k) +
                                            kappa[0] * (r + force) -
                                            kappa[0] * mu::d_dx_P(model, i, j, k);
                    ;
                }
            }
        }

        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceV(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_V(model, nu, i, j, k);
                    rhs_buff.V(i, j, k) = (r + force);

                    model_buff.V(i, j, k) = model.V(i, j, k) +
                                            kappa[0] * (r + force) -
                                            kappa[0] * mu::d_dy_P(model, i, j, k);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceW(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_W(model, nu, i, j, k);
                    rhs_buff.W(i, j, k) = (r + force);

                    model_buff.W(i, j, k) = model.W(i, j, k) +
                                            kappa[0] * (r + force) -
                                            kappa[0] * mu::d_dz_P(model, i, j, k);
                }
            }
        }

        boundary_cond.apply(model_buff, time + kappa[0]);
    }

    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        // COMPUTE B TERM
        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {
                    b[i + (j + k * ny) * nx] = kappa[7] * (mu::d_dx_U(model_buff, i, j, k) +
                                                           mu::d_dy_V(model_buff, i, j, k) +
                                                           mu::d_dz_W(model_buff, i, j, k));
                }
            }
        }
        // solver.setb(b);
        // SOLVE FOR phi2-pn
        // solver.solve(X);
        // UPDATE PRESSURE (not real pressure, but phi2-pn, used to compute the gradient)
        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {
                    // model_buff.P(i, j, k) = X[i + (j + k * ny) * nx];
                }
            }
        }
    }

    /// Y2 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        //         for (index_t k = 0; k < nz; ++k)
        //         {
        //             for (index_t j = 0; j < ny; ++j)
        //             {
        // #pragma omp simd
        //                 for (index_t i = 0; i < nx; ++i)
        //                 {
        //                     model.U(i, j, k) = model_buff.U(i, j, k) /*- mu::d_dx_P(model_buff, i, j, k) / kappa[7]*/;
        //                     model.V(i, j, k) = model_buff.V(i, j, k) /*- mu::d_dy_P(model_buff, i, j, k) / kappa[7]*/;
        //                     model.W(i, j, k) = model_buff.W(i, j, k) /*- mu::d_dz_P(model_buff, i, j, k) / kappa[7]*/;
        //                     // Also update pressure with the real value, i.e. phi2 (model holds pn, model_buff holds phi2-pn)
        //                     // model.P(i, j, k) += model_buff.P(i, j, k) ;
        //                 }
        //             }
        //         }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef ForcingT
    ft.set_time(time + kappa[0]);
#endif

    /// Y3* //////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force2 = getForceU(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r1 = rhs_buff.U(i, j, k);
                    const Real r2 = rhs_U(model_buff, nu, i, j, k);

                    rhs_buff.U(i, j, k) = (r2 + force2);

                    model.U(i, j, k) = model_buff.U(i, j, k) -
                                       kappa[1] * r1 +
                                       kappa[2] * (r2 + force2) -
                                       kappa[5] * mu::d_dx_P(model_buff, i, j, k);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force2 = getForceV(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r1 = rhs_buff.V(i, j, k);
                    const Real r2 = rhs_V(model_buff, nu, i, j, k);

                    rhs_buff.V(i, j, k) = (r2 + force2);

                    model.V(i, j, k) = model_buff.V(i, j, k) -
                                       kappa[1] * r1 +
                                       kappa[2] * (r2 + force2) -
                                       kappa[5] * mu::d_dy_P(model_buff, i, j, k);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force2 = getForceW(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r1 = rhs_buff.W(i, j, k);
                    const Real r2 = rhs_W(model_buff, nu, i, j, k);

                    rhs_buff.W(i, j, k) = (r2 + force2);

                    model.W(i, j, k) = model_buff.W(i, j, k) -
                                       kappa[1] * r1 +
                                       kappa[2] * (r2 + force2) -
                                       kappa[5] * mu::d_dz_P(model_buff, i, j, k);
                }
            }
        }

        boundary_cond.apply(model, time + kappa[4]);
    }

    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        // COMPUTE B TERM
        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {
                    b[i + (j + k * ny) * nx] = kappa[8] * (mu::d_dx_U(model_buff, i, j, k) +
                                                           mu::d_dy_V(model_buff, i, j, k) +
                                                           mu::d_dz_W(model_buff, i, j, k));
                }
            }
        }
        // solver.setb(b);
        // SOLVE FOR phi3-phi2
        // solver.solve(X);
        // UPDATE PRESSURE (not real pressure, but phi3-phi2, used to compute the gradient)
        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {
                    // model_buff.P(i, j, k) = X[i + (j + k * ny) * nx];
                }
            }
        }
    }
    /// Y3 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        //         for (index_t k = 0; k < nz; ++k)
        //         {
        //             for (index_t j = 0; j < ny; ++j)
        //             {
        // #pragma omp simd
        //                 for (index_t i = 0; i < nx; ++i)
        //                 {
        //                     model.U(i, j, k) = model_buff.U(i, j, k) /*- mu::d_dx_P(model_buff, i, j, k) / kappa[8]*/;
        //                     model.V(i, j, k) = model_buff.V(i, j, k) /*- mu::d_dy_P(model_buff, i, j, k) / kappa[8]*/;
        //                     model.W(i, j, k) = model_buff.W(i, j, k) /*- mu::d_dz_P(model_buff, i, j, k) / kappa[8]*/;
        //                     // Also update pressure with the real value, i.e. phi3 (model holds phi2, model_buff holds phi3-phi2)
        //                     // model.P(i, j, k) += model_buff.P(i, j, k);
        //                 }
        //             }
        //         }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef ForcingT
    ft.set_time(time + kappa[4]);
#endif

    /// u(n+1)* //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceU(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_U(model, nu, i, j, k);

                    model_buff.U(i, j, k) = model.U(i, j, k) -
                                            kappa[2] * rhs_buff.U(i, j, k) +
                                            kappa[3] * (r + force) -
                                            kappa[6] * mu::d_dx_P(model, i, j, k);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceV(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_V(model, nu, i, j, k);

                    model_buff.V(i, j, k) = model.V(i, j, k) -
                                            kappa[2] * rhs_buff.V(i, j, k) +
                                            kappa[3] * (r + force) -
                                            kappa[6] * mu::d_dy_P(model, i, j, k);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceW(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_W(model, nu, i, j, k);

                    model_buff.W(i, j, k) = model.W(i, j, k) -
                                            kappa[2] * rhs_buff.W(i, j, k) +
                                            kappa[3] * (r + force) -
                                            kappa[6] * mu::d_dz_P(model, i, j, k);
                }
            }
        }

        boundary_cond.apply(model_buff, time + deltat);
    }

    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        // COMPUTE B TERM
        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {
                    b[i + (j + k * ny) * nx] = kappa[9] * (mu::d_dx_U(model_buff, i, j, k) +
                                                           mu::d_dy_V(model_buff, i, j, k) +
                                                           mu::d_dz_W(model_buff, i, j, k));
                }
            }
        }
        // solver.setb(b);
        // SOLVE FOR pn+1-phi3
        // solver.solve(X);
        // UPDATE PRESSURE (not real pressure, but pn+1-phi3, used to compute the gradient)
        for (index_t k = 0; k < nz; ++k)
        {
            for (index_t j = 0; j < ny; ++j)
            {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i)
                {
                    // model_buff.P(i, j, k) = X[i + (j + k * ny) * nx];
                }
            }
        }
    }

    /// un+1 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        //         for (index_t k = 0; k < nz; ++k)
        //         {
        //             for (index_t j = 0; j < ny; ++j)
        //             {
        // #pragma omp simd
        //                 for (index_t i = 0; i < nx; ++i)
        //                 {
        //                     model.U(i, j, k) = model_buff.U(i, j, k) /*- mu::d_dx_P(model_buff, i, j, k) / kappa[9]*/;
        //                     model.V(i, j, k) = model_buff.V(i, j, k) /*- mu::d_dy_P(model_buff, i, j, k) / kappa[9]*/;
        //                     model.W(i, j, k) = model_buff.W(i, j, k) /*- mu::d_dz_P(model_buff, i, j, k) / kappa[9]*/;
        //                     // Also update pressure with the real value, i.e. phi3 (model holds phi2, model_buff holds phi3-phi2)
        //                     // model.P(i, j, k) += model_buff.P(i, j, k);
        //                 }
        //             }
        //         }
    }
    model.swap(model_buff);
}

#endif // AEROHPC_A_RUNGEKUTTA_H
