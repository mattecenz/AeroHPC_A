#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include "GridData.hpp"
#include "ForcingTerm.hpp"
#include "Boundaries.hpp"
#include "mathUtils.hpp"
#include "poissonSolver.hpp"


///TODO TO REMOVE ONLY FOR DEBUG
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
///////////////////////////////

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
void rungeKutta(GridData &model, GridData &model_buff,
                GridData &rhs_buff, GridData &pressure_buff,
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

    // Convert pressure_buff from velocity convention to X and b convention
#define X_at(i,j,k) pressure_buff.U(i,j,k)
#define b_at(i,j,k) pressure_buff.V(i,j,k)
#define X (&pressure_buff.U(0,0,0))
#define b (&pressure_buff.V(0,0,0))

    /// Y2* //////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceU(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_U(model, nu, i, j, k);
                    rhs_buff.U(i, j, k) = (r + force);

                    model_buff.U(i, j, k) = model.U(i, j, k)
                                            + k_0 * (r + force)
#ifdef ENABLE_PRESSURE
                                            - k_0 * mu::dp_dx_U(model, i, j, k)
#endif
                            ;
                }
            }
        }

        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
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

                    model_buff.V(i, j, k) = model.V(i, j, k)
                                            + k_0 * (r + force)
#ifdef ENABLE_PRESSURE
                                            - k_0 * mu::dp_dy_V(model, i, j, k)
#endif
                            ;
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
                    rhs_buff.W(i, j, k) = (r + force);

                    model_buff.W(i, j, k) = model.W(i, j, k)
                                            + k_0 * (r + force)
#ifdef ENABLE_PRESSURE
                                            - k_0 * mu::dp_dz_W(model, i, j, k)
#endif
                            ;
                }
            }
        }

        b_print(model_buff, dir, 1);

        boundary_cond.apply(model_buff, t_0);

        b_print(model_buff, dir, 2);
    }

#ifdef ENABLE_PRESSURE
    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        // COMPUTE B TERM
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    b_at(i, j, k) = inv_k_0 * mu::vel_div(model_buff, i, j, k);
                }
            }
        }

        p_solver.setB(b);

        // SOLVE FOR phi2-pn
        p_solver.solve(X);
    }

    /// UPDATE PRESSURE /////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
                memcpy(&model_buff.P(0, j, k), &X_at(0, j, k), sizeof(Real) * nx);
            }
        }

        b_print(model_buff, dir, 3);

        boundary_cond.apply(model_buff, t_0);

        b_print(model_buff, dir, 4);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /// Y2 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    model_buff.U(i, j, k) -= k_0 * mu::dp_dx_U(model_buff, i, j, k);
                    model_buff.V(i, j, k) -= k_0 * mu::dp_dy_V(model_buff, i, j, k);
                    model_buff.W(i, j, k) -= k_0 * mu::dp_dz_W(model_buff, i, j, k);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    model_buff.P(i, j, k) += model.P(i, j, k);
                }
            }
        }

        b_print(model_buff, dir, 5);

        boundary_cond.apply(model_buff, t_0);

        b_print(model_buff, dir, 6);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif

#ifdef ForcingT
    ft.set_time(t_0);
#endif

    /// Y3* //////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
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
                                       - k_1 * r1
                                       + k_2 * (r2 + force2)
#ifdef ENABLE_PRESSURE
                                       - k_3 * mu::dp_dx_U(model_buff, i, j, k)
#endif
                            ;
                }
            }
        }

        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
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
                                       - k_1 * r1
                                       + k_2 * (r2 + force2)
#ifdef ENABLE_PRESSURE
                                       - k_3 * mu::dp_dy_V(model_buff, i, j, k)
#endif
                            ;
                }
            }
        }

        for (index_t k = 0; k < nz; ++k) {
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
                                       - k_1 * r1
                                       + k_2 * (r2 + force2)
#ifdef ENABLE_PRESSURE
                                       - k_3 * mu::dp_dz_W(model_buff, i, j, k)
#endif
                            ;
                }
            }
        }

        b_print(model, dir, 7);

        boundary_cond.apply(model, t_1);

        b_print(model, dir, 8);
    }

#ifdef ENABLE_PRESSURE
    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        // COMPUTE B TERM
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    b_at(i, j, k) = inv_k_3 * mu::vel_div(model, i, j, k);
                }
            }
        }

        p_solver.setB(b);

        // SOLVE FOR phi3-phi2
        p_solver.solve(X);
    }

    /// UPDATE PRESSURE /////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
                memcpy(&model.P(0, j, k), &X_at(0, j, k), sizeof(Real) * nx);
            }
        }

        b_print(model, dir, 9);

        boundary_cond.apply(model, t_1);

        b_print(model, dir, 10);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Y3 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    model.U(i, j, k) -= k_3 * mu::dp_dx_U(model, i, j, k);
                    model.V(i, j, k) -= k_3 * mu::dp_dy_V(model, i, j, k);
                    model.W(i, j, k) -= k_3 * mu::dp_dz_W(model, i, j, k);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    model.P(i, j, k) += model_buff.P(i, j, k);
                }
            }
        }

        b_print(model, dir, 11);

        boundary_cond.apply(model, t_1);

        b_print(model, dir, 12);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif


#ifdef ForcingT
    ft.set_time(t_1);
#endif

    /// u(n+1)* //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
#ifdef ForcingT
                    getPhys(i, j, k);
                    const Real force = getForceU(ft);
#else
                    constexpr Real force = 0;
#endif

                    const Real r = rhs_U(model, nu, i, j, k);

                    model_buff.U(i, j, k) = model.U(i, j, k)
                                            - k_4 * rhs_buff.U(i, j, k)
                                            + k_5 * (r + force)
#ifdef ENABLE_PRESSURE
                                            - k_6 * mu::dp_dx_U(model, i, j, k)
#endif
                            ;
                }
            }
        }

        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
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
                                            - k_4 * rhs_buff.V(i, j, k)
                                            + k_5 * (r + force)
#ifdef ENABLE_PRESSURE
                                            - k_6 * mu::dp_dy_V(model, i, j, k)
#endif
                            ;
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
                                            - k_4 * rhs_buff.W(i, j, k)
                                            + k_5 * (r + force)
#ifdef ENABLE_PRESSURE
                                            - k_6 * mu::dp_dz_W(model, i, j, k)
#endif
                            ;
                }
            }
        }

        b_print(model_buff, dir, 13);

        boundary_cond.apply(model_buff, t_2);

        b_print(model_buff, dir, 14);
    }

#ifdef ENABLE_PRESSURE
    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        // COMPUTE B TERM
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    b_at(i, j, k) = inv_k_6 * mu::vel_div(model_buff, i, j, k);
                }
            }
        }

        p_solver.setB(b);


        // SOLVE FOR pn+1-phi3
        p_solver.solve(X);
    }

    /// UPDATE PRESSURE /////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
                memcpy(&model_buff.P(0, j, k), &X_at(0, j, k), sizeof(Real) * nx);
            }
        }

        b_print(model_buff, dir, 15);

        boundary_cond.apply(model_buff, t_2);

        b_print(model_buff, dir, 16);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// un+1 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    model_buff.U(i, j, k) -= k_6 * mu::dp_dx_U(model_buff, i, j, k);
                    model_buff.V(i, j, k) -= k_6 * mu::dp_dy_V(model_buff, i, j, k);
                    model_buff.W(i, j, k) -= k_6 * mu::dp_dz_W(model_buff, i, j, k);
                }
            }
        }

        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    model_buff.P(i, j, k) += model.P(i, j, k);
                }
            }
        }

        b_print(model_buff, dir, 17);

        boundary_cond.apply(model_buff, t_2);

        b_print(model_buff, dir, 18);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif

    model.swap(model_buff);
}

#endif // AEROHPC_A_RUNGEKUTTA_H
