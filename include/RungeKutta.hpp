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

#define getForce(ftt) {  ftt.computeGx(px + dx, py + sdy, pz + sdz),     \
                        ftt.computeGy(px + sdx, py + dy, pz + sdz),      \
                        ftt.computeGz(px + sdx, py + sdy, pz + dz)   }
#endif

using namespace utils;


struct RKConst {
    static constexpr Real alpha0 = 8.0 / 15.0;
    static constexpr Real alpha1 = 17.0 / 60.0;
    static constexpr Real alpha2 = 5.0 / 12.0;
    static constexpr Real alpha3 = 3.0 / 4.0;
    static constexpr Real alpha4 = 2.0 / 3.0;
};

//RHS function
template<Addressing_T A>
inline Vector rhs(const Grid<A> &grid,const Real nu,const index_t i,const index_t j,const index_t k) {
    return - conv(grid,i,j,k) + nu * lap(grid,i,j,k);
}

template<Component C>
inline Real rhss(Grid<STANDARD> &grid, const Real nu, const index_t i,const  index_t j,const  index_t k) {
    return - conv<C>(grid,i,j,k) + nu * lap<C>(grid,i,j,k);
}


inline Real rhss_u(Grid<STANDARD> &grid, const Real nu, const index_t i,const  index_t j,const  index_t k) {
    return - conv_u(grid,i,j,k) + nu * lap_u(grid,i,j,k);
}




//Runge-Kutta method
void rungeKutta(Grid<STANDARD> &model, Grid<STANDARD> &model_buff, Grid<STANDARD> &rhs_buff,
                Real reynolds, Real deltat, Real time,
                Boundaries<STANDARD> &boundary_cond) {

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

    //grid -> Y1
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

    //Y2.

    for (index_t i = 0; i < nx; ++i) {
        for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
            for (index_t k = 0; k < nz; ++k) {

                const Real r = rhss_u(model, nu, i, j, k);

#ifdef ForcingT
                getPhys(i, j, k);

                Real force = ft.computeGx(px + dx, py + sdy, pz + sdz);

                rhs_buff(U,i,j,k) = (r + force);

                //TODO BRANCHED CODE IS THE DEVIL
                //if (i != nx - 1)
                    model_buff(U, i, j, k) = model(U, i, j, k) + kappa[0] * (r + force);
#else
                Y2(U, i, j, k) = model(U, i, j, k) + kappa[0] * r[0];
                Y2(V, i, j, k) = model(V, i, j, k) + kappa[0] * r[1];
                Y2(W, i, j, k) = model(W, i, j, k) + kappa[0] * r[2];
#endif

            }

        }
    }
    for (index_t i = 0; i < nx; ++i) {
        for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
            for (index_t k = 0; k < nz; ++k) {

                const Real r = rhss<V>(model, nu, i, j, k);

#ifdef ForcingT
                getPhys(i, j, k);

                Real force = ft.computeGy(px + sdx, py + dy, pz + sdz);

                rhs_buff(V,i,j,k) = (r + force);

                //TODO BRANCHED CODE IS THE DEVIL
                //if (i != nx - 1)
                model_buff(V, i, j, k) = model(V, i, j, k) + kappa[0] * (r + force);
#else
                Y2(U, i, j, k) = model(U, i, j, k) + kappa[0] * r[0];
                Y2(V, i, j, k) = model(V, i, j, k) + kappa[0] * r[1];
                Y2(W, i, j, k) = model(W, i, j, k) + kappa[0] * r[2];
#endif

            }

        }
    }

    for (index_t i = 0; i < nx; ++i) {
        for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
            for (index_t k = 0; k < nz; ++k) {

                const Real r = rhss<W>(model, nu, i, j, k);

#ifdef ForcingT
                getPhys(i, j, k);

                Real force = ft.computeGz(px + sdx, py + sdy, pz + dz);

                rhs_buff(W,i,j,k) = (r + force);

                //TODO BRANCHED CODE IS THE DEVIL
                //if (i != nx - 1)
                model_buff(W, i, j, k) = model(W, i, j, k) + kappa[0] * (r + force);
#else
                Y2(U, i, j, k) = model(U, i, j, k) + kappa[0] * r[0];
                Y2(V, i, j, k) = model(V, i, j, k) + kappa[0] * r[1];
                Y2(W, i, j, k) = model(W, i, j, k) + kappa[0] * r[2];
#endif

            }

        }
    }


    boundary_cond.apply(model_buff, time + kappa[0]);

#ifdef ForcingT
     ft.set_time(time + kappa[0]);
#endif

    //Y3.

    for (index_t i = 0; i < nx; ++i) {
        for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
            for (index_t k = 0; k < nz; ++k) {

                Vector r2 = rhs(model_buff, nu, i, j, k);

#ifdef ForcingT

                getPhys(i, j, k);

                const Vector force2 = getForce(ft);

                const Real rhs_u = rhs_buff(U, i, j, k);
                const Real rhs_v = rhs_buff(V, i, j, k);
                const Real rhs_w = rhs_buff(W, i, j, k);

                rhs_buff(U,i,j,k) = (r2[0] + force2[0]);
                rhs_buff(V,i,j,k) = (r2[1] + force2[1]);
                rhs_buff(W,i,j,k) = (r2[2] + force2[2]);

                // TODO BRANCHED CODE IS THE DEVIL
                //if (i != nx - 1)
                    model(U, i, j, k) = model_buff(U, i, j, k)
                                          - kappa[1] * rhs_u
                                          + kappa[2] * (r2[0] + force2[0]);
                //if (j != ny - 1)
                    model(V, i, j, k) = model_buff(V, i, j, k) +
                                          -kappa[1] * rhs_v
                                          + kappa[2] * (r2[1] + force2[1]);
                //if (k != nz - 1)
                    model(W, i, j, k) = model_buff(W, i, j, k) +
                                          -kappa[1] * rhs_w
                                          + kappa[2] * (r2[2] + force2[2]);
#else
                Y3(U, i, j, k) = Y2(U, i, j, k)
                                        - kappa[1] * r[0]
                                        + kappa[2] * r2[0];
                Y3(V, i, j, k) = Y2(V, i, j, k) +
                                        - kappa[1] * r[1]
                                        + kappa[2] * r2[1];
                Y3(W, i, j, k) = Y2(W, i, j, k) +
                                        - kappa[1] * r[2]
                                        + kappa[2] * r2[2];
#endif


            }
        }
    }


    boundary_cond.apply(model, time + kappa[4]);

#ifdef ForcingT
    ft.set_time(time + kappa[4]);
#endif

    //u(n+1)

    for (index_t i = 0; i < nx; ++i) {
        for (index_t j = 0; j < ny; ++j) {
#pragma omp simd
            for (index_t k = 0; k < nz; ++k) {
                const Vector r = rhs(model, nu, i, j, k);

#ifdef ForcingT
                getPhys(i, j, k);

                // Be careful because fts are now inverted
                Vector force = getForce(ft);

                // TODO BRANCHED CODE IS THE DEVIL
                //if (i != nx - 1)
                    model_buff(U, i, j, k) = model(U, i, j, k)
                                             - kappa[2] * rhs_buff(U,i,j,k)
                                             + kappa[3] * (r[0] + force[0]);
                //if (j != ny - 1)
                    model_buff(V, i, j, k) = model(V, i, j, k) +
                                             -kappa[2] * rhs_buff(V,i,j,k)
                                             + kappa[3] * (r[1] + force[1]);
                //if (k != nz - 1)
                    model_buff(W, i, j, k) = model(W, i, j, k) +
                                             -kappa[2] * rhs_buff(W,i,j,k)
                                             + kappa[3] * (r[2] + force[2]);
#else
                model(U, i, j, k) = Y3(U, i, j, k)
                                        - kappa[2] * r[0]
                                        + kappa[3] * r2[0];
                model(V, i, j, k) = Y3(V, i, j, k) +
                                        - kappa[2] * r[1]
                                        + kappa[3] * r2[1];
                model(W, i, j, k) = Y3(W, i, j, k) +
                                        - kappa[2] * r[2]
                                        + kappa[3] * r2[2];
#endif


            }
        }
    }

    boundary_cond.apply(model_buff, time + deltat);

    model.swap(model_buff);


}


#endif //AEROHPC_A_RUNGEKUTTA_H
