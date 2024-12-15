#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include "GridData.hpp"
#include "ForcingTerm.hpp"
#include "Boundaries.hpp"
#include "mathUtils.hpp"
// #include "PoissonSolver.hpp"

#ifdef ForcingT
#define getPhys(i, j, k)  Real px = real(i + model.structure.px) * dx;   \
                        Real py = real(j + model.structure.py) * dy;     \
                        Real pz = real(k + model.structure.pz) * dz

#define getForceU(ftt) ftt.computeGx(px + dx, py + sdy, pz + sdz)

#define getForceV(ftt) ftt.computeGy(px + sdx, py + dy, pz + sdz)

#define getForceW(ftt) ftt.computeGz(px + sdx, py + sdy, pz + dz)
#endif

namespace mu = mathUtils;

struct RKConst {
    // Constants for velocity terms
    static constexpr Real alpha0 = 64.0 / 120.0;
    static constexpr Real alpha1 = 34.0 / 120.0;
    static constexpr Real alpha2 = 50.0 / 120.0;
    static constexpr Real alpha3 = 90.0 / 120.0;
    static constexpr Real alpha4 = 2.0 / 3.0;
    // Constants for pressure terms
    static constexpr Real beta0 = 64.0 / 120.0;
    static constexpr Real beta1 = 16.0 / 120.0;
    static constexpr Real beta2 = 40.0 / 120.0;
    // Constants for FastPoissonSolver terms
    static constexpr Real gamma0 = 120.0 / 64.0;
    static constexpr Real gamma1 = 64.0 / 120.0;
    static constexpr Real gamma2 = 40.0 / 120.0;
};

////RHS function
#define rhs(C) inline Real rhs_##C(GridData &grid, const Real nu, const index_t i, const index_t j, const index_t k) { \
    return -mu::conv_##C(grid, i, j, k) + nu * mu::lap_##C(grid, i, j, k); \
}

rhs(U)

rhs(V)

rhs(W)


//Runge-Kutta method
void rungeKutta(GridData &model, GridData &model_buff, GridData &rhs_buff,
                Real reynolds, Real deltat, Real time,
                Boundaries &boundary_cond, 
                C2Decomp *c2d) {

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

    // bring this out later
    double L = 1.0;
    // Define the solver object
    // poissonSolver solver(nx, ny, nz, L, c2d);

    //kappa -> weighted_deltat 
    std::array<const Real, 11> kappa{
            // Constants for velocity terms
            RKConst::alpha0 * deltat,
            RKConst::alpha1 * deltat,
            RKConst::alpha2 * deltat,
            RKConst::alpha3 * deltat,
            RKConst::alpha4 * deltat,
            // Constants for pressure terms
            RKConst::beta0 * deltat,
            RKConst::beta1 * deltat,
            RKConst::beta2 * deltat,
            // Constants for FastPoissonSolver terms
            RKConst::gamma0 / deltat,
            RKConst::gamma1 * deltat,
            RKConst::gamma2 * deltat
    };

#ifdef ForcingT
    ForcingTerm ft(reynolds, time);
#endif

/// Y2 //////////////////////////////////////////////////////////////////////////////////////////////
    {
    // First Stage - Prediction Step: Calculating Y2* (Velocity Prediction)

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

                    model_buff.U(i, j, k) = model.U(i, j, k) + kappa[0] * (r + force) - kappa[5]*mathUtils::d_dx_P(model, i, j, k);
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

                    model_buff.V(i, j, k) = model.V(i, j, k) + kappa[0] * (r + force) - kappa[5] * mathUtils::d_dy_P(model, i, j, k);
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

                    model_buff.W(i, j, k) = model.W(i, j, k) + kappa[0] * (r + force) - kappa[5] * mathUtils::d_dz_P(model, i, j, k);
                }
            }
        }


        boundary_cond.apply(model_buff, time + kappa[0]);
    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// COLLECT DATA ///////////////////////////////////////////////////////////////////////////////////////////////
//TODO: THIS CAN BE COMPUTED IN THE RUNGE-KUTTA ITERATIONS, BUT LOOP ORDERS HAVE TO BE CHANGED
Real b[nx * ny * nz];

for (index_t k = 0; k < nz; ++k) {
    for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
        for (index_t i = 0; i < nx; ++i) {
            b[i + (j + k * ny) * nx] =  kappa[8] * 
                                            (
                                                mathUtils::d_dx_U(model_buff,i,j,k)*model_buff.U(i,j,k)+
                                                mathUtils::d_dy_V(model_buff,i,j,k)*model_buff.V(i,j,k)+
                                                mathUtils::d_dz_W(model_buff,i,j,k)*model_buff.W(i,j,k)
                                            );
        }
    }
}

Real x[nx * ny * nz];
// TODO:
    // CALL SOLVER
    // COMPUTE NEW PRESSURE AND PUT IT IN THE GRID
    // ADJUST VELOCITY

    // solver.setB(b);
    // solver.solve(x);
    
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
                                       + kappa[2] * (r2 + force2)
                                       - kappa[6] * mathUtils::d_dx_P(model_buff,i,j,k);
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
                                       + kappa[2] * (r2 + force2)
                                       - kappa[6] * mathUtils::d_dy_P(model_buff,i,j,k);
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
                                       + kappa[2] * (r2 + force2)
                                       - kappa[6] * mathUtils::d_dz_P(model_buff,i,j,k);
                }
            }
        }

        boundary_cond.apply(model, time + kappa[4]);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// COLLECT DATA ///////////////////////////////////////////////////////////////////////////////////////////////
//TODO: THIS CAN BE COMPUTED IN THE RUNGE-KUTTA ITERATIONS, BUT LOOP ORDERS HAVE TO BE CHANGED
for (index_t k = 0; k < nz; ++k) {
    for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
        for (index_t i = 0; i < nx; ++i) {
            b[i + (j + k * ny) * nx] =  kappa[8] * 
                                            (
                                                mathUtils::d_dx_U(model,i,j,k)*model.U(i,j,k)+
                                                mathUtils::d_dy_V(model,i,j,k)*model.V(i,j,k)+
                                                mathUtils::d_dz_W(model,i,j,k)*model.W(i,j,k)
                                            );
        }
    }
}
//TODO:
    // CALL SOLVER
    // COMPUTE NEW PRESSURE AND PUT IT IN THE GRID
    // ADJUST VELOCITY
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
                                            + kappa[3] * (r + force)
                                            - kappa[7] * mathUtils::d_dx_P(model,i,j,k);
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
                                            - kappa[7] * mathUtils::d_dy_P(model,i,j,k);
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
                                            - kappa[7] * mathUtils::d_dz_P(model,i,j,k);
                }
            }
        }

        boundary_cond.apply(model_buff, time + deltat);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// COLLECT DATA ///////////////////////////////////////////////////////////////////////////////////////////////
//TODO: THIS CAN BE COMPUTED IN THE RUNGE-KUTTA ITERATIONS, BUT LOOP ORDERS HAVE TO BE CHANGED
for (index_t k = 0; k < nz; ++k) {
    for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
        for (index_t i = 0; i < nx; ++i) {
            b[i + (j + k * ny) * nx] =  kappa[8] * 
                                            (
                                                mathUtils::d_dx_U(model_buff,i,j,k)*model_buff.U(i,j,k)+
                                                mathUtils::d_dy_V(model_buff,i,j,k)*model_buff.V(i,j,k)+
                                                mathUtils::d_dz_W(model_buff,i,j,k)*model_buff.W(i,j,k)
                                            );
        }
    }
}

// TODO:
    // CALL SOLVER
    // COMPUTE NEW PRESSURE AND PUT IT IN THE GRID
    // ADJUST VELOCITY
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    model.swap(model_buff);
}


#endif //AEROHPC_A_RUNGEKUTTA_H