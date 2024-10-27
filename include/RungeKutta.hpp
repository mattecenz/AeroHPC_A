#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include <StaggeredGrid.hpp>
#include <ForcingTerm.hpp>
#include <chrono>
#include <utils.hpp>
#include <Model.hpp>

#define getPhys(i, j, k)  Real px = real(i) * model.dx;\
                        Real py = real(j) * model.dy;\
                        Real pz = real(k) * model.dz\


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
inline Vector rhs(Model<A> &model, index_t i, index_t j, index_t k) {

    Real nu = (1 / model.reynolds);

    Vector convs = conv(model, i, j, k);

    Vector laps{
            lap(model, U, i, j, k),
            lap(model, V, i, j, k),
            lap(model, W, i, j, k)
    };

    return -convs + nu * laps;
}

//Runge-Kutta method
void rungeKutta(Model<STANDARD> &model, Model<STANDARD> &Y2, Model<STANDARD> &Y3, Real deltat, Real time) {

    Real sdx = model.sdx;
    Real sdy = model.sdy;
    Real sdz = model.sdz;

    //grid -> Y1
    //kappa -> weighted_deltat 
    std::array<Real, 5> kappa{
            RKConst::alpha0 * deltat,
            RKConst::alpha1 * deltat,
            RKConst::alpha2 * deltat,
            RKConst::alpha3 * deltat,
            RKConst::alpha4 * deltat
    };

#ifdef ForcingT
    ForcingTerm ft(model.reynolds, time);
#endif

    // Then here changes as our indexes should be changed a bit
    // But careful, on the face i=0 we need to solve only the x velocity
    // Its better if we duplicate some code or else we will go crazy
    for (index_t j = 1; j < model.grid.ny-1; ++j) {
        for (index_t k = 1; k < model.grid.nz-1; ++k) {
            Vector r = rhs(model, 0, j, k);

#ifdef ForcingT
            getPhys(0, j, k);

            Vector force{
                    ft.computeGx(px + sdx, py, pz), //Gx
                    ft.computeGy(px, py + sdy, pz), //Gy
                    ft.computeGz(px, py, pz + sdz)  //Gz
            };

            Y2.grid(U, 0, j, k) = model.grid(U, 0, j, k) + kappa[0] * (r[0] + force[0]);
#else
            Y2.grid(U, 0, j, k) = model.grid(U, 0, j, k) + kappa[0] * r[0];
#endif
        }
    }
    // And its ugly but we need to do the same for the two faces also at y=0 and z=0
    // Maybe wrap in some function or think how to optimize
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t k = 1; k < model.grid.nz-1; ++k) {
            Vector r = rhs(model, i, 0, k);

#ifdef ForcingT
            getPhys(i, 0, k);

            Vector force{
                    ft.computeGx(px + sdx, py, pz), //Gx
                    ft.computeGy(px, py + sdy, pz), //Gy
                    ft.computeGz(px, py, pz + sdz)  //Gz
            };

            Y2.grid(V, i, 0, k) = model.grid(V, i, 0, k) + kappa[0] * (r[1] + force[1]);
#else
            Y2.grid(V, i, 0, k) = model.grid(V, i, 0, k) + kappa[0] * r[1];
#endif
        }
    }
    // z = 0
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t j = 1; j < model.grid.ny-1; ++j) {
            Vector r = rhs(model, i, j, 0);

#ifdef ForcingT
            getPhys(i, j, 0);

            Vector force{
                    ft.computeGx(px + sdx, py, pz), //Gx
                    ft.computeGy(px, py + sdy, pz), //Gy
                    ft.computeGz(px, py, pz + sdz)  //Gz
            };

            Y2.grid(W, i, j, 0) = model.grid(W, i, j, 0) + kappa[0] * (r[2] + force[2]);
#else
            Y2.grid(W, i, j, 0) = model.grid(W, i, j, 0) + kappa[0] * r[2];
#endif
        }
    }
    //Y2.
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t j = 1; j < model.grid.ny-1; ++j) {
            for (index_t k = 1; k < model.grid.nz-1; ++k) {
                Vector r = rhs(model, i, j, k);

#ifdef ForcingT
                getPhys(i, j, k);

                Vector force{
                    ft.computeGx(px + sdx, py, pz), //Gx
                    ft.computeGy(px, py + sdy, pz), //Gy
                    ft.computeGz(px, py, pz + sdz)  //Gz
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

    Y2.applyBCs(time + kappa[0]);

#ifdef ForcingT
    ForcingTerm ft2(model.reynolds, time + kappa[0]);
#endif

    // We kinda need to do the same ugly stuff also for Y3
    // x = 0
    for (index_t j = 1; j < model.grid.ny-1; ++j) {
        for (index_t k = 1; k < model.grid.nz-1; ++k) {
            Vector r = rhs(model, 0, j, k);
            Vector r2 = rhs(Y2, 0, j, k);

#ifdef ForcingT
            getPhys(0, j, k);

            Vector force{
                    ft.computeGx(px + sdx, py, pz), //Gx
                    ft.computeGy(px, py + sdy, pz), //Gy
                    ft.computeGz(px, py, pz + sdz)  //Gz
            };

            Vector force2{
                    ft2.computeGx(px + sdx, py, pz), //Gx
                    ft2.computeGy(px, py + sdy, pz), //Gy
                    ft2.computeGz(px, py, pz + sdz)  //Gz
            };

            Y3.grid(U, 0, j, k) = Y2.grid(U, 0, j, k)
                                        - kappa[1] * (r[0] + force[0])
                                        + kappa[2] * (r2[0] + force2[0]);
#else
            Y3.grid(U, 0, j, k) = Y2.grid(U, 0, j, k)
                                        - kappa[1] * r[0]
                                        + kappa[2] * r2[0];
#endif
        }
    }
    // y = 0
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t k = 1; k < model.grid.nz-1; ++k) {
            Vector r = rhs(model, i, 0, k);
            Vector r2 = rhs(Y2, i, 0, k);

#ifdef ForcingT
            getPhys(i, 0, k);

            Vector force{
                    ft.computeGx(px + sdx, py, pz), //Gx
                    ft.computeGy(px, py + sdy, pz), //Gy
                    ft.computeGz(px, py, pz + sdz)  //Gz
            };

            Vector force2{
                    ft2.computeGx(px + sdx, py, pz), //Gx
                    ft2.computeGy(px, py + sdy, pz), //Gy
                    ft2.computeGz(px, py, pz + sdz)  //Gz
            };

            Y3.grid(V, i, 0, k) = Y2.grid(V, i, 0, k) +
                                        -kappa[1] * (r[1] + force[1])
                                        + kappa[2] * (r2[1] + force2[1]);
#else
            Y3.grid(V, i, 0, k) = Y2.grid(V, i, 0, k)
                                    - kappa[1] * r[0]
                                    + kappa[2] * r2[0];
#endif
        }
    }
    // z = 0
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t j = 1; j < model.grid.ny-1; ++j) {
            Vector r = rhs(model, i, j, 0);
            Vector r2 = rhs(Y2, i, j, 0);

#ifdef ForcingT
            getPhys(i, j, 0);

            Vector force{
                    ft.computeGx(px + sdx, py, pz), //Gx
                    ft.computeGy(px, py + sdy, pz), //Gy
                    ft.computeGz(px, py, pz + sdz)  //Gz
            };

            Vector force2{
                    ft2.computeGx(px + sdx, py, pz), //Gx
                    ft2.computeGy(px, py + sdy, pz), //Gy
                    ft2.computeGz(px, py, pz + sdz)  //Gz
            };

            Y3.grid(W, i, j, 0) = Y2.grid(W, i, j, 0) +
                                        -kappa[1] * (r[2] + force[2])
                                        + kappa[2] * (r2[2] + force2[2]);
#else
            Y3.grid(W, i, j, 0) = Y2.grid(W, i, j, 0) +
                                    - kappa[1] * r[2]
                                    + kappa[2] * r2[2];
#endif
        }
    }
    //Y3.
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t j = 1; j < model.grid.ny-1; ++j) {
            for (index_t k = 1; k < model.grid.nz-1; ++k) {
                Vector r = rhs(model, i, j, k);
                Vector r2 = rhs(Y2, i, j, k);

#ifdef ForcingT

                getPhys(i, j, k);

                Vector force{
                    ft.computeGx(px + sdx, py, pz), //Gx
                    ft.computeGy(px, py + sdy, pz), //Gy
                    ft.computeGz(px, py, pz + sdz)  //Gz
            };

            Vector force2{
                    ft2.computeGx(px + sdx, py, pz), //Gx
                    ft2.computeGy(px, py + sdy, pz), //Gy
                    ft2.computeGz(px, py, pz + sdz)  //Gz
            };

                Y3.grid(U, i, j, k) = Y2.grid(U, i, j, k)
                                        - kappa[1] * (r[0] + force[0])
                                        + kappa[2] * (r2[0] + force2[0]);
                Y3.grid(V, i, j, k) = Y2.grid(V, i, j, k) +
                                        -kappa[1] * (r[1] + force[1])
                                        + kappa[2] * (r2[1] + force2[1]);
                Y3.grid(W, i, j, k) = Y2.grid(W, i, j, k) +
                                        -kappa[1] * (r[2] + force[2])
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


    Y3.applyBCs(time + kappa[4]);

#ifdef ForcingT
    ft.set_time(time + kappa[4]);
#endif

    // And surprise surprise we need to do also the same shit here
    // x = 0
    for (index_t j = 1; j < model.grid.ny-1; ++j) {
        for (index_t k = 1; k < model.grid.nz-1; ++k) {
            Vector r = rhs(Y2, 0, j, k);
            Vector r2 = rhs(Y3, 0, j, k);

#ifdef ForcingT
            getPhys(0, j, k);

            Vector force{
                    ft2.computeGx(px + sdx, py, pz), //Gx
                    ft2.computeGy(px, py + sdy, pz), //Gy
                    ft2.computeGz(px, py, pz + sdz)  //Gz
            };

            Vector force2{
                    ft.computeGx(px + sdx, py, pz), //Gx
                    ft.computeGy(px, py + sdy, pz), //Gy
                    ft.computeGz(px, py, pz + sdz)  //Gz
            };

            model.grid(U, 0, j, k) = Y3.grid(U, 0, j, k)
                                            - kappa[2] * (r[0] + force[0])
                                            + kappa[3] * (r2[0] + force2[0]);
#else
            model.grid(U, 0, j, k) = Y3.grid(U, 0, j, k)
                                    - kappa[2] * r[0]
                                    + kappa[3] * r2[0];
#endif
        }
    }
    // y = 0
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t k = 1; k < model.grid.nz-1; ++k) {
            Vector r = rhs(Y2, i, 0, k);
            Vector r2 = rhs(Y3, i, 0, k);

#ifdef ForcingT
            getPhys(i, 0, k);

            Vector force{
                    ft2.computeGx(px + sdx, py, pz), //Gx
                    ft2.computeGy(px, py + sdy, pz), //Gy
                    ft2.computeGz(px, py, pz + sdz)  //Gz
            };

            Vector force2{
                    ft.computeGx(px + sdx, py, pz), //Gx
                    ft.computeGy(px, py + sdy, pz), //Gy
                    ft.computeGz(px, py, pz + sdz)  //Gz
            };

            model.grid(V, i, 0, k) = Y3.grid(V, i, 0, k) +
                                        -kappa[2] * (r[1] + force[1])
                                        + kappa[3] * (r2[1] + force2[1]);
#else
            model.grid(V, i, 0, k) = Y3.grid(V, i, 0, k) +
                                    - kappa[2] * r[1]
                                    + kappa[3] * r2[1];
#endif
        }
    }
    // z = 0
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t j = 1; j < model.grid.ny-1; ++j) {
            Vector r = rhs(Y2, i, j, 0);
            Vector r2 = rhs(Y3, i, j, 0);

#ifdef ForcingT
            getPhys(i, j, 0);

            Vector force{
                    ft2.computeGx(px + sdx, py, pz), //Gx
                    ft2.computeGy(px, py + sdy, pz), //Gy
                    ft2.computeGz(px, py, pz + sdz)  //Gz
            };

            Vector force2{
                    ft.computeGx(px + sdx, py, pz), //Gx
                    ft.computeGy(px, py + sdy, pz), //Gy
                    ft.computeGz(px, py, pz + sdz)  //Gz
            };

            model.grid(W, i, j, 0) = Y3.grid(W, i, j, 0) +
                                        -kappa[2] * (r[2] + force[2])
                                        + kappa[3] * (r2[2] + force2[2]);
#else
            model.grid(W, i, j, 0) = Y3.grid(W, i, j, 0) +
                                    - kappa[2] * r[2]
                                    + kappa[3] * r2[2];
#endif
        }
    }
    //u(n+1)    
    for (index_t i = 1; i < model.grid.nx-1; ++i) {
        for (index_t j = 1; j < model.grid.ny-1; ++j) {
            for (index_t k = 1; k < model.grid.nz-1; ++k) {
                Vector r = rhs(Y2, i, j, k);
                Vector r2 = rhs(Y3, i, j, k);

#ifdef ForcingT
                getPhys(i, j, k);

                // Be careful because fts are now inverted
                Vector force{
                        ft2.computeGx(px + sdx, py, pz), //Gx
                        ft2.computeGy(px, py + sdy, pz), //Gy
                        ft2.computeGz(px, py, pz + sdz)  //Gz
                };

                Vector force2{
                        ft.computeGx(px + sdx, py, pz), //Gx
                        ft.computeGy(px, py + sdy, pz), //Gy
                        ft.computeGz(px, py, pz + sdz)  //Gz
                };

                model.grid(U, i, j, k) = Y3.grid(U, i, j, k)
                                            - kappa[2] * (r[0] + force[0])
                                            + kappa[3] * (r2[0] + force2[0]);
                model.grid(V, i, j, k) = Y3.grid(V, i, j, k) +
                                            -kappa[2] * (r[1] + force[1])
                                            + kappa[3] * (r2[1] + force2[1]);
                model.grid(W, i, j, k) = Y3.grid(W, i, j, k) +
                                            -kappa[2] * (r[2] + force[2])
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

    model.applyBCs(time + deltat);

}


#endif //AEROHPC_A_RUNGEKUTTA_H
