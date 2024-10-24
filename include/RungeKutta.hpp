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
                int i, int j, int k,
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
void rungeKutta(int n, Real Re, Model<Addressing_T::STANDARD> model,
                Model<Addressing_T::STANDARD> &model_out, double deltat, double time) {

    //grid -> Y1
    //kappa -> weighted_deltat 
    std::array<double, 5> kappa{
            RKConstants::alpha1 * deltat,
            RKConstants::alpha2 * deltat,
            RKConstants::alpha3 * deltat,
            RKConstants::alpha4 * deltat,
            RKConstants::alpha5 * deltat
    };

    //BUFFERS
    Model<Addressing_T::STANDARD> Y2(model);
    Model<Addressing_T::STANDARD> Y3(model);

#ifdef ForcingT
    ForcingTerm forcingterm(Re);
#endif

    //Y2.
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {

#ifdef ForcingT
                Y2.grid(Component::U, i, j, k) = model.grid(Component::U, i, j, k) +
                                                 rhs(model, Component::U, i, j, k, kappa[0], Re, forcingterm, time);
                Y2.grid(Component::V, i, j, k) = model.grid(Component::V, i, j, k) +
                                                 rhs(model, Component::V, i, j, k, kappa[0], Re, forcingterm, time);
                Y2.grid(Component::W, i, j, k) = model.grid(Component::W, i, j, k) +
                                                 rhs(model, Component::W, i, j, k, kappa[0], Re, forcingterm, time);
#else
                Y2.grid(Component::U, i, j, k) = model.grid(Component::U, i, j, k) + rhs(model, Component::U, i, j, k, kappa[0], Re);
                Y2.grid(Component::V, i, j, k) = model.grid(Component::V, i, j, k) + rhs(model, Component::V, i, j, k, kappa[0], Re);
                Y2.grid(Component::W, i, j, k) = model.grid(Component::W, i, j, k) + rhs(model, Component::W, i, j, k, kappa[0], Re);
#endif

            }
        }
    }

    //Y3.
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {

#ifdef ForcingT
                Y3.grid(Component::U, i, j, k) = Y2.grid(Component::U, i, j, k) +
                                                 rhs(Y2, Component::U, i, j, k, kappa[1], Re, forcingterm,
                                                     (time + kappa[0])) -
                                                 rhs(model, Component::U, i, j, k, kappa[2], Re, forcingterm, time);
                Y3.grid(Component::V, i, j, k) = Y2.grid(Component::V, i, j, k) +
                                                 rhs(Y2, Component::V, i, j, k, kappa[1], Re, forcingterm,
                                                     (time + kappa[0])) -
                                                 rhs(model, Component::V, i, j, k, kappa[2], Re, forcingterm, time);
                Y3.grid(Component::W, i, j, k) = Y2.grid(Component::W, i, j, k) +
                                                 rhs(Y2, Component::W, i, j, k, kappa[1], Re, forcingterm,
                                                     (time + kappa[0])) -
                                                 rhs(model, Component::W, i, j, k, kappa[2], Re, forcingterm, time);
#else
                Y3.grid(Component::U, i, j, k) = Y2.grid(Component::U, i, j, k) + rhs(Y2, Component::U, i, j, k, kappa[1], Re) -
                                            rhs(model, Component::U, i, j, k, kappa[2], Re);
                Y3.grid(Component::V, i, j, k) = Y2.grid(Component::V, i, j, k) + rhs(Y2, Component::V, i, j, k, kappa[1], Re) -
                                            rhs(model, Component::V, i, j, k, kappa[2], Re);
                Y3.grid(Component::W, i, j, k) = Y2.grid(Component::W, i, j, k) + rhs(Y2, Component::W, i, j, k, kappa[1], Re) -
                                            rhs(model, Component::W, i, j, k, kappa[2], Re);
#endif


            }
        }
    }

    //u(n+1)    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {

#ifdef ForcingT
                model_out.grid(Component::U, i, j, k) =
                        Y3.grid(Component::U, i, j, k) +
                        rhs(Y3, Component::U, i, j, k, kappa[3], Re, forcingterm, (time + kappa[4])) -
                        rhs(Y2, Component::U, i, j, k, kappa[1], Re, forcingterm, (time + kappa[0]));
                model_out.grid(Component::V, i, j, k) =
                        Y3.grid(Component::V, i, j, k) +
                        rhs(Y3, Component::V, i, j, k, kappa[3], Re, forcingterm, (time + kappa[4])) -
                        rhs(Y2, Component::V, i, j, k, kappa[1], Re, forcingterm, (time + kappa[0]));
                model_out.grid(Component::W, i, j, k) =
                        Y3.grid(Component::W, i, j, k) +
                        rhs(Y3, Component::W, i, j, k, kappa[3], Re, forcingterm, (time + kappa[4])) -
                        rhs(Y2, Component::W, i, j, k, kappa[1], Re, forcingterm, (time + kappa[0]));
#else
                model_out.grid(Component::U, i, j, k) =
                        Y3.grid(Component::U, i, j, k) + rhs(Y3, Component::U, i, j, k, kappa[3], Re) -
                        rhs(Y2, Component::U, i, j, k, kappa[1], Re);
                model_out.grid(Component::V, i, j, k) =
                        Y3.grid(Component::V, i, j, k) + rhs(Y3, Component::V, i, j, k, kappa[3], Re) -
                        rhs(Y2, Component::V, i, j, k, kappa[1], Re);
                model_out.grid(Component::W, i, j, k) =
                        Y3.grid(Component::W, i, j, k) + rhs(Y3, Component::W, i, j, k, kappa[3], Re) -
                        rhs(Y2, Component::W, i, j, k, kappa[1], Re);
#endif


            }
        }
    }

}


#endif //AEROHPC_A_RUNGEKUTTA_H
