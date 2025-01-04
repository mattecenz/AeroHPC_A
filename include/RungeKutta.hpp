#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include "GridData.hpp"
#include "ForcingTerm.hpp"
#include "Boundaries.hpp"
#include "mathUtils.hpp"
#include "poissonSolver.hpp"

#include "printBuffer.hpp"


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
    static constexpr Real alpha0 = 8.0 / 15.0;
    static constexpr Real alpha1 = 17.0 / 60.0;
    static constexpr Real alpha2 = 5.0 / 12.0;
    static constexpr Real alpha3 = 3.0 / 4.0;
    static constexpr Real alpha4 = 2.0 / 3.0;
    static constexpr Real alpha5 = 2.0 / 15.0;
    static constexpr Real alpha6 = 1.0 / 3.0;
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
    if (!rank) {
        create_directories(dir);
    }

    if (!rank) {
        std::string nn{dir + "0"};
        print(model, nn);
    }

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
    std::array<const Real, 8> kappa{
        RKConst::alpha0 * deltat,
        RKConst::alpha1 * deltat,
        RKConst::alpha2 * deltat,
        RKConst::alpha3 * deltat,
        RKConst::alpha4 * deltat,
        RKConst::alpha5 * deltat,
        RKConst::alpha6 * deltat,
        real(1.0 / (RKConst::alpha6 * deltat))
    };

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
                                            + kappa[0] * (r + force)
                                            - kappa[0] * mu::d_dx_P(model, i, j, k);
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
                                            + kappa[0] * (r + force)
                                            - kappa[0] * mu::d_dy_P(model, i, j, k);
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
                                            + kappa[0] * (r + force)
                                            - kappa[0] * mu::d_dz_P(model, i, j, k);
                }
            }
        }

        if (!rank) {
            std::string nn{dir + "1"};
            print(model_buff, nn);
        }

        boundary_cond.apply(model_buff, time + kappa[0]);

        if (!rank) {
            std::string nn{dir + "2"};
            print(model_buff, nn);
        }
    }

    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        // COMPUTE B TERM
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    b_at(i, j, k) = kappa[7] * (mu::d_dx_U(model_buff, i, j, k)
                                                + mu::d_dy_V(model_buff, i, j, k)
                                                + mu::d_dz_W(model_buff, i, j, k));
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

        if (!rank) {
            std::string nn{dir + "3"};
            print(model_buff, nn);
        }

        boundary_cond.apply(model_buff, time + kappa[0]);

        if (!rank) {
            std::string nn{dir + "4"};
            print(model_buff, nn);
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /// Y2 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    model_buff.U(i, j, k) -= kappa[0] * mu::d_dx_P(model_buff, i, j, k);
                    model_buff.V(i, j, k) -= kappa[0] * mu::d_dy_P(model_buff, i, j, k);
                    model_buff.W(i, j, k) -= kappa[0] * mu::d_dz_P(model_buff, i, j, k);

                    model_buff.P(i, j, k) += model.P(i, j, k);
                }
            }
        }

        if (!rank) {
            std::string nn{dir + "5"};
            print(model_buff, nn);
        }

        boundary_cond.apply(model_buff, time + kappa[0]);

        if (!rank) {
            std::string nn{dir + "6"};
            print(model_buff, nn);
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef ForcingT
    ft.set_time(time + kappa[0]);
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
                                       - kappa[1] * r1
                                       + kappa[2] * (r2 + force2)
                                       - kappa[5] * mu::d_dx_P(model_buff, i, j, k);
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
                                       - kappa[1] * r1
                                       + kappa[2] * (r2 + force2)
                                       - kappa[5] * mu::d_dy_P(model_buff, i, j, k);
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
                                       - kappa[1] * r1
                                       + kappa[2] * (r2 + force2)
                                       - kappa[5] * mu::d_dz_P(model_buff, i, j, k);
                }
            }
        }

        if (!rank) {
            std::string nn{dir + "7"};
            print(model, nn);
        }
        boundary_cond.apply(model, time + kappa[4]);

        if (!rank) {
            std::string nn{dir + "8"};
            print(model, nn);
        }
    }

    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        // COMPUTE B TERM
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    b_at(i, j, k) = kappa[7] * (mu::d_dx_U(model, i, j, k)
                                                + mu::d_dy_V(model, i, j, k)
                                                + mu::d_dz_W(model, i, j, k));
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

        if (!rank) {
            std::string nn{dir + "9"};
            print(model, nn);
        }

        boundary_cond.apply(model, time + kappa[4]);

        if (!rank) {
            std::string nn{dir + "10"};
            print(model, nn);
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Y3 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    model.U(i, j, k) -= kappa[0] * mu::d_dx_P(model, i, j, k);
                    model.V(i, j, k) -= kappa[0] * mu::d_dy_P(model, i, j, k);
                    model.W(i, j, k) -= kappa[0] * mu::d_dz_P(model, i, j, k);
                    model.P(i, j, k) += model_buff.P(i, j, k);
                }
            }
        }

        if (!rank) {
            std::string nn{dir + "11"};
            print(model, nn);
        }

        boundary_cond.apply(model, time + kappa[4]);

        if (!rank) {
            std::string nn{dir + "12"};
            print(model, nn);
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef ForcingT
    ft.set_time(time + kappa[4]);
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
                                            - kappa[2] * rhs_buff.U(i, j, k) +
                                            + kappa[3] * (r + force)
                                            - kappa[6] * mu::d_dx_P(model, i, j, k);
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
                                            - kappa[2] * rhs_buff.V(i, j, k)
                                            + kappa[3] * (r + force)
                                            - kappa[6] * mu::d_dy_P(model, i, j, k);
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
                                            + kappa[3] * (r + force)
                                            - kappa[6] * mu::d_dz_P(model, i, j, k);
                }
            }
        }

        if (!rank) {
            std::string nn{dir + "13"};
            print(model_buff, nn);
        }

        boundary_cond.apply(model_buff, time + deltat);

        if (!rank) {
            std::string nn{dir + "14"};
            print(model_buff, nn);
        }
    }

    /// POISSON SOLVER ///////////////////////////////////////////////////////////////////////////////////
    {
        // COMPUTE B TERM
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    b_at(i, j, k) = kappa[7] * (mu::d_dx_U(model_buff, i, j, k) +
                                                mu::d_dy_V(model_buff, i, j, k) +
                                                mu::d_dz_W(model_buff, i, j, k));
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

        if (!rank) {
            std::string nn{dir + "15"};
            print(model_buff, nn);
        }

        boundary_cond.apply(model_buff, time + deltat);

        if (!rank) {
            std::string nn{dir + "16"};
            print(model_buff, nn);
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// un+1 /////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        for (index_t k = 0; k < nz; ++k) {
            for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
                for (index_t i = 0; i < nx; ++i) {
                    model_buff.U(i, j, k) -= kappa[6] * mu::d_dx_P(model_buff, i, j, k);
                    model_buff.V(i, j, k) -= kappa[6] * mu::d_dy_P(model_buff, i, j, k);
                    model_buff.W(i, j, k) -= kappa[6] * mu::d_dz_P(model_buff, i, j, k);
                    model_buff.P(i, j, k) += model.P(i, j, k);
                }
            }
        }

        if (!rank) {
            std::string nn{dir + "17"};
            print(model_buff, nn);
        }

        boundary_cond.apply(model_buff, time + deltat);

        if (!rank) {
            std::string nn{dir + "18"};
            print(model_buff, nn);
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    model.swap(model_buff);
}

#endif // AEROHPC_A_RUNGEKUTTA_H
