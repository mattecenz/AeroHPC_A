#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include <StaggeredGrid.hpp>
#include <ForcingTerm.hpp>
#include <chrono>
#include <utils.hpp>
#include <Model.hpp>

using namespace utils;


struct RKConstants {
    static constexpr double alpha1 = 64.0 / 120.0;
    static constexpr double alpha2 = 50.0 / 120.0;
    static constexpr double alpha3 = 34.0 / 120.0;
    static constexpr double alpha4 = 90.0 / 120.0;
    static constexpr double alpha5 = 2.0 / 3.0;
};

//RHS function
template<Addressing_T A>
inline Real rhs(Model<A> &model, Component c,
                index_t i, index_t j, index_t k,
                double Kappa, Real Re
#ifdef ForcingT
                , ForcingTerm &forcingterm, double time
#endif
) {

#ifdef ForcingT
    forcingterm.update_time(time);
#endif

    return Kappa * (
            -(conv(model, c, i, j, k)) + (1 / Re) * lap(model, c, i, j, k)
            #ifdef ForcingT
            + forcingterm.compute(i, j, k, c)
#endif
    );
}

//Runge-Kutta method
void rungeKutta(Real Re, Model<STANDARD> &model, Model<STANDARD> &Y2, Model<STANDARD> &Y3, double deltat, double time) {

    //grid -> Y1
    //kappa -> weighted_deltat 
    std::array<double, 5> kappa{
            RKConstants::alpha1 * deltat,
            RKConstants::alpha2 * deltat,
            RKConstants::alpha3 * deltat,
            RKConstants::alpha4 * deltat,
            RKConstants::alpha5 * deltat
    };

#ifdef ForcingT
    ForcingTerm forcingterm(Re);
#endif

    //Y2.
    for (index_t i = 0; i < model.grid.nx; ++i) {
        for (index_t j = 0; j < model.grid.ny; ++j) {
            for (index_t k = 0; k < model.grid.nz; ++k) {

#ifdef ForcingT
                Y2.grid(U, i, j, k) = model.grid(U, i, j, k) +
                                                 rhs(model, U, i, j, k, kappa[0], Re, forcingterm, time);
                Y2.grid(V, i, j, k) = model.grid(V, i, j, k) +
                                                 rhs(model, V, i, j, k, kappa[0], Re, forcingterm, time);
                Y2.grid(W, i, j, k) = model.grid(W, i, j, k) +
                                                 rhs(model, W, i, j, k, kappa[0], Re, forcingterm, time);
#else
                Y2.grid(U, i, j, k) = model.grid(U, i, j, k) + rhs(model, U, i, j, k, kappa[0], Re);
                Y2.grid(V, i, j, k) = model.grid(V, i, j, k) + rhs(model, V, i, j, k, kappa[0], Re);
                Y2.grid(W, i, j, k) = model.grid(W, i, j, k) + rhs(model, W, i, j, k, kappa[0], Re);
#endif

            }
        }
    }

    //Y3.
    for (index_t i = 0; i < model.grid.nx; ++i) {
        for (index_t j = 0; j < model.grid.ny; ++j) {
            for (index_t k = 0; k < model.grid.nz; ++k) {

#ifdef ForcingT
                Y3.grid(U, i, j, k) = Y2.grid(U, i, j, k) +
                                                 rhs(Y2, U, i, j, k, kappa[1], Re, forcingterm,
                                                     (time + kappa[0])) -
                                                 rhs(model, U, i, j, k, kappa[2], Re, forcingterm, time);
                Y3.grid(V, i, j, k) = Y2.grid(V, i, j, k) +
                                                 rhs(Y2, V, i, j, k, kappa[1], Re, forcingterm,
                                                     (time + kappa[0])) -
                                                 rhs(model, V, i, j, k, kappa[2], Re, forcingterm, time);
                Y3.grid(W, i, j, k) = Y2.grid(W, i, j, k) +
                                                 rhs(Y2, W, i, j, k, kappa[1], Re, forcingterm,
                                                     (time + kappa[0])) -
                                                 rhs(model, W, i, j, k, kappa[2], Re, forcingterm, time);
#else
                Y3.grid(U, i, j, k) = Y2.grid(U, i, j, k) + rhs(Y2, U, i, j, k, kappa[1], Re) -
                                            rhs(model, U, i, j, k, kappa[2], Re);
                Y3.grid(V, i, j, k) = Y2.grid(V, i, j, k) + rhs(Y2, V, i, j, k, kappa[1], Re) -
                                            rhs(model, V, i, j, k, kappa[2], Re);
                Y3.grid(W, i, j, k) = Y2.grid(W, i, j, k) + rhs(Y2, W, i, j, k, kappa[1], Re) -
                                            rhs(model, W, i, j, k, kappa[2], Re);
#endif


            }
        }
    }

    //u(n+1)    
    for (index_t i = 0; i < model.grid.nx; ++i) {
        for (index_t j = 0; j < model.grid.ny; ++j) {
            for (index_t k = 0; k < model.grid.nz; ++k) {

#ifdef ForcingT
                model.grid(U, i, j, k) =
                        Y3.grid(U, i, j, k) +
                        rhs(Y3, U, i, j, k, kappa[3], Re, forcingterm, (time + kappa[4])) -
                        rhs(Y2, U, i, j, k, kappa[1], Re, forcingterm, (time + kappa[0]));
                model.grid(V, i, j, k) =
                        Y3.grid(V, i, j, k) +
                        rhs(Y3, V, i, j, k, kappa[3], Re, forcingterm, (time + kappa[4])) -
                        rhs(Y2, V, i, j, k, kappa[1], Re, forcingterm, (time + kappa[0]));
                model.grid(W, i, j, k) =
                        Y3.grid(W, i, j, k) +
                        rhs(Y3, W, i, j, k, kappa[3], Re, forcingterm, (time + kappa[4])) -
                        rhs(Y2, W, i, j, k, kappa[1], Re, forcingterm, (time + kappa[0]));
#else
                model_out.grid(U, i, j, k) =
                        Y3.grid(U, i, j, k) + rhs(Y3, U, i, j, k, kappa[3], Re) -
                        rhs(Y2, U, i, j, k, kappa[1], Re);
                model_out.grid(V, i, j, k) =
                        Y3.grid(V, i, j, k) + rhs(Y3, V, i, j, k, kappa[3], Re) -
                        rhs(Y2, V, i, j, k, kappa[1], Re);
                model_out.grid(W, i, j, k) =
                        Y3.grid(W, i, j, k) + rhs(Y3, W, i, j, k, kappa[3], Re) -
                        rhs(Y2, W, i, j, k, kappa[1], Re);
#endif


            }
        }
    }

}


#endif //AEROHPC_A_RUNGEKUTTA_H
