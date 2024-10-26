#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include <StaggeredGrid.hpp>
#include <ForcingTerm.hpp>
#include <chrono>
#include <utils.hpp>
#include <Model.hpp>

using namespace utils;


struct RKConst {
    static constexpr Real alpha0 = 8.0 / 15.0;
    static constexpr Real alpha1 = 17.0 / 60.0;
    static constexpr Real alpha2 = 5.0 / 12.0;
    static constexpr Real alpha3 = 3.0 / 4.0;
#ifdef ForcingT
    static constexpr Real alpha4 = 2.0 / 3.0;
#endif
};

//RHS function
template<Addressing_T A>
inline Vector rhs(Model<A> &model, index_t i, index_t j, index_t k) {

    Real nu = (1 / model.reynolds);

    Vector convs = conv(model, i, j, k);

    Vector laps{
            lap(model, U, i, j, k),
            lap(model, V, i, j, k),
            lap(model, W, i, j, k)
    };

    return - convs + nu * laps;
}

//Runge-Kutta method
void rungeKutta(Model<STANDARD> &model, Model<STANDARD> &Y2, Model<STANDARD> &Y3, Real deltat, Real time) {

    Real sdx = model.sdx;
    Real sdy = model.sdy;
    Real sdz = model.sdz;

    //grid -> Y1
    //kappa -> weighted_deltat 
    std::array<Real, 4> kappa{
            RKConst::alpha0 * deltat,
            RKConst::alpha1 * deltat,
            RKConst::alpha2 * deltat,
            RKConst::alpha3 * deltat,
    };

#ifdef ForcingT
    Real kappa4 = RKConst::alpha4 * deltat;

    ForcingTerm ft(model.reynolds, time);
#endif

    //Y2.
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t j = 1; j < model.grid.ny-1; ++j) {
            for (index_t k = 1; k < model.grid.nz-1; ++k) {
                Vector r = rhs(model, i, j, k);

#ifdef ForcingT
                Real px = static_cast<Real>(i) * model.dx;
                Real py = static_cast<Real>(j) * model.dy;
                Real pz = static_cast<Real>(k) * model.dz;

                Vector force{
                        ft.computeGx(px + sdx, py, pz), //Gx
                        ft.computeGy(px, py + sdy, pz), //Gy
                        ft.computeGz(px, py, pz + sdz) //Gz
                };

                Y2.grid(U, i, j, k) = model.grid(U, i, j, k) + kappa[0] * (r[0] + force[0]);
                Y2.grid(V, i, j, k) = model.grid(V, i, j, k) + kappa[0] * (r[1] + force[1]);
                Y2.grid(W, i, j, k) = model.grid(W, i, j, k) + kappa[0] * (r[2] + force[2]);
#else
                Y2.grid(U, i, j, k) = model.grid(U, i, j, k) + kappa[0] * r[0];
                Y2.grid(V, i, j, k) = model.grid(V, i, j, k) + kappa[0] * r[1];
                Y2.grid(W, i, j, k) = model.grid(W, i, j, k) + kappa[0] * r[2];
#endif

            }
        }
    }

#ifdef ForcingT
    ForcingTerm ft2(model.reynolds, time + kappa[0]);
#endif

    //Y3.
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t j = 1; j < model.grid.ny-1; ++j) {
            for (index_t k = 1; k < model.grid.nz-1; ++k) {
                Vector r = rhs(model, i, j, k);
                Vector r2 = rhs(Y2, i, j, k);

#ifdef ForcingT

                Real px = static_cast<Real>(i) * model.dx;
                Real py = static_cast<Real>(j) * model.dy;
                Real pz = static_cast<Real>(k) * model.dz;

                Vector force{
                        ft.computeGx(px + sdx, py, pz),
                        ft.computeGy(px, py + sdy, pz),
                        ft.computeGz(px, py, pz + sdz)
                };

                Vector force2{
                        ft2.computeGx(px + sdx, py, pz),
                        ft2.computeGy(px, py + sdy, pz),
                        ft2.computeGz(px, py, pz + sdz)
                };


                Y3.grid(U, i, j, k) = Y2.grid(U, i, j, k)
                                        - kappa[1] * (r[0] + force[0])
                                        + kappa[2] * (r2[0] + force2[0]);
                Y3.grid(V, i, j, k) = Y2.grid(V, i, j, k) +
                                        - kappa[1] * (r[1] + force[1])
                                        + kappa[2] * (r2[1] + force2[1]);
                Y3.grid(W, i, j, k) = Y2.grid(W, i, j, k) +
                                        - kappa[1] * (r[2] + force[2])
                                        + kappa[2] * (r2[2] + force2[2]);
#else
                Y3.grid(U, i, j, k) = Y2.grid(U, i, j, k)
                                        - kappa[1] * r[0]
                                        + kappa[2] * r2[0];
                Y3.grid(V, i, j, k) = Y2.grid(V, i, j, k) +
                                        - kappa[1] * r[1]
                                        + kappa[2] * r2[1];
                Y3.grid(W, i, j, k) = Y2.grid(W, i, j, k) +
                                        - kappa[1] * r[2]
                                        + kappa[2] * r2[2];
#endif


            }
        }
    }


#ifdef ForcingT
    ft.set_time(time + kappa4);
#endif

    //u(n+1)    
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t j = 1; j < model.grid.ny-1; ++j) {
            for (index_t k = 1; k < model.grid.nz-1; ++k) {
                Vector r = rhs(Y2, i, j, k);
                Vector r2 = rhs(Y3, i, j, k);

#ifdef ForcingT
                Real px = static_cast<Real>(i) * model.dx;
                Real py = static_cast<Real>(j) * model.dy;
                Real pz = static_cast<Real>(k) * model.dz;

                // Be careful because fts are now inverted
                Vector force{
                        ft2.computeGx(px + sdx, py, pz),
                        ft2.computeGy(px, py + sdy, pz),
                        ft2.computeGz(px, py, pz + sdz)
                };

                Vector force2{
                        ft.computeGx(px + sdx, py, pz),
                        ft.computeGy(px, py + sdy, pz),
                        ft.computeGz(px, py, pz + sdz)
                };


                model.grid(U, i, j, k) = Y3.grid(U, i, j, k)
                                      - kappa[2] * (r[0] + force[0])
                                      + kappa[3] * (r2[0] + force2[0]);
                model.grid(V, i, j, k) = Y3.grid(V, i, j, k) +
                                      - kappa[2] * (r[1] + force[1])
                                      + kappa[3] * (r2[1] + force2[1]);
                model.grid(W, i, j, k) = Y3.grid(W, i, j, k) +
                                      - kappa[2] * (r[2] + force[2])
                                      + kappa[3] * (r2[2] + force2[2]);
#else
                model.grid(U, i, j, k) = Y3.grid(U, i, j, k)
                                        - kappa[2] * r[0]
                                        + kappa[3] * r2[0];
                model.grid(V, i, j, k) = Y3.grid(V, i, j, k) +
                                        - kappa[2] * r[1]
                                        + kappa[3] * r2[1];
                model.grid(W, i, j, k) = Y3.grid(W, i, j, k) +
                                        - kappa[2] * r[2]
                                        + kappa[3] * r2[2];
#endif


            }
        }
    }

}


#endif //AEROHPC_A_RUNGEKUTTA_H
