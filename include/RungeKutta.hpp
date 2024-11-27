#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include "Grid.hpp"
#include "ForcingTerm.hpp"
#include "Boundaries.hpp"
#include "utils.hpp"

#ifdef ForcingT
#define getPhys(i, j, k)  Real px = real(i) * dx;   \
                        Real py = real(j) * dy;     \
                        Real pz = real(k) * dz

#define getForceU(ftt) ftt.computeGx(px + dx, py + sdy, pz + sdz)

#define getForceV(ftt) ftt.computeGy(px + sdx, py + dy, pz + sdz)

#define getForceW(ftt) ftt.computeGz(px + sdx, py + sdy, pz + dz)

#endif

using namespace utils;


struct RKConst {
    static constexpr Real alpha0 = 8.0 / 15.0;
    static constexpr Real alpha1 = 17.0 / 60.0;
    static constexpr Real alpha2 = 5.0 / 12.0;
    static constexpr Real alpha3 = 3.0 / 4.0;
    static constexpr Real alpha4 = 2.0 / 3.0;
};

////RHS function
#define rhs(C) inline Real rhs_##C(Grid &grid, const Real nu, const index_t i, const index_t j, const index_t k) { \
    return -conv_##C(grid, i, j, k) + nu * lap_##C(grid, i, j, k); \
}

rhs(U)

rhs(V)

rhs(W)


//Runge-Kutta method
void rungeKutta(Grid &model, Grid &model_buff, Grid &rhs_buff,
                Real reynolds, Real deltat, Real time,
                Boundaries &boundary_cond) {

    const Real nu = (real(1) / reynolds);

    const Real dx = model.dx;
    const Real dy = model.dy;
    const Real dz = model.dz;

    const Real sdx = model.sdx;
    const Real sdy = model.sdy;
    const Real sdz = model.sdz;

    const index_t nx = model.nx;
    const index_t ny = model.ny;
    const index_t nz = model.nz;

    //kappa -> weighted_deltat 
    std::array<const Real, 5> kappa{
            RKConst::alpha0 * deltat,
            RKConst::alpha1 * deltat,
            RKConst::alpha2 * deltat,
            RKConst::alpha3 * deltat,
            RKConst::alpha4 * deltat
    };

#ifdef ForcingT
    ForcingTerm ft(reynolds, time);
#endif

/// Y2 //////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx - 1; ++i) {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceU(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_U(model, nu, i, j, k);
                    rhs_buff.U(i, j, k) = (r + force);

                    model_buff.U(i, j, k) = model.U(i, j, k) + kappa[0] * (r + force);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny - 1; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceV(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_V(model, nu, i, j, k);
                    rhs_buff.V(i, j, k) = (r + force);

                    model_buff.V(i, j, k) = model.V(i, j, k) + kappa[0] * (r + force);
                }

            }
        }

        for (index_t k = 0; k < nz - 1; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceW(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_W(model, nu, i, j, k);
                    rhs_buff.W(i, j, k) = (r + force);

                    model_buff.W(i, j, k) = model.W(i, j, k) + kappa[0] * (r + force);
                }
            }
        }


        boundary_cond.apply(model_buff, time + kappa[0]);
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef ForcingT
    ft.set_time(time + kappa[0]);
#endif

/// Y3 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx - 1; ++i) {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force2 = getForceU(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r1 = rhs_buff.U(i, j, k);
                    const Real r2 = rhs_U(model_buff, nu, i, j, k);

                    rhs_buff.U(i, j, k) = (r2 + force2);

                    model.U(i, j, k) = model_buff.U(i, j, k)
                                       - kappa[1] * r1
                                       + kappa[2] * (r2 + force2);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny - 1; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force2 = getForceV(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r1 = rhs_buff.V(i, j, k);
                    const Real r2 = rhs_V(model_buff, nu, i, j, k);

                    rhs_buff.V(i, j, k) = (r2 + force2);

                    model.V(i, j, k) = model_buff.V(i, j, k)
                                       - kappa[1] * r1
                                       + kappa[2] * (r2 + force2);
                }
            }
        }

        for (index_t k = 0; k < nz - 1; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force2 = getForceW(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r1 = rhs_buff.W(i, j, k);
                    const Real r2 = rhs_W(model_buff, nu, i, j, k);

                    rhs_buff.W(i, j, k) = (r2 + force2);

                    model.W(i, j, k) = model_buff.W(i, j, k)
                                       - kappa[1] * r1
                                       + kappa[2] * (r2 + force2);
                }
            }
        }

        boundary_cond.apply(model, time + kappa[4]);
    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef ForcingT
    ft.set_time(time + kappa[4]);
#endif

/// u(n+1) //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx - 1; ++i) {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceU(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_U(model, nu, i, j, k);

                    model_buff.U(i, j, k) = model.U(i, j, k)
                                            - kappa[2] * rhs_buff.U(i, j, k)
                                            + kappa[3] * (r + force);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny - 1; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceV(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_V(model, nu, i, j, k);

                    model_buff.V(i, j, k) = model.V(i, j, k)
                                            - kappa[2] * rhs_buff.V(i, j, k)
                                            + kappa[3] * (r + force);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {

#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceW(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_W(model, nu, i, j, k);

                    model_buff.W(i, j, k) = model.W(i, j, k)
                                            - kappa[2] * rhs_buff.W(i, j, k)
                                            + kappa[3] * (r + force);
                }
            }
        }

        boundary_cond.apply(model_buff, time + deltat);
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    model.swap(model_buff);


}


#endif //AEROHPC_A_RUNGEKUTTA_H
