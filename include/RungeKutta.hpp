#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include <StaggeredGrid.hpp>
#include <ForcingTerm.hpp>
#include <chrono>
#include <utils.hpp>

using namespace utils;


struct RKConstants {
    static constexpr double alpha1 = 64.0 / 120.0;
    static constexpr double alpha2 = 50.0 / 120.0;
    static constexpr double alpha3 = 34.0 / 120.0;
    static constexpr double alpha4 = 90.0 / 120.0;
    static constexpr double alpha5 = 2.0 / 3.0;
};

#ifdef ForcingTerm

inline Real rhsforcingterm(const ForcingTerm &forcingterm, Component c, int i, int j, int k, double timeUpdate){
    forcingterm.update_time(timeUpdate);
    return forcingterm.compute(i,j,k,c);
}

//RHS function
template<Addressing_T A>
inline Real rhs(const StaggeredGrid<A> &grid, Component c, int i, int j, int k, double Kappa, Real Re  
,const ForcingTerm &forcingterm, double time) {
    return Kappa * (-(conv(grid, c, i, j, k)) + (1 / Re) * lap(grid, c, i, j, k)    
    + rhsforcingterm(forcingterm,c,i,j,k,time));
}

#endif

//RHS function
template<Addressing_T A>
inline Real rhs(const StaggeredGrid<A> &grid, Component c, int i, int j, int k, double Kappa, Real Re) {
    return Kappa * (-(conv(grid, c, i, j, k)) + (1 / Re) * lap(grid, c, i, j, k));
}


//Runge-Kutta method
void rungeKutta(int n, Real Re, const StaggeredGrid<Addressing_T::STANDARD> grid,
                StaggeredGrid<Addressing_T::STANDARD> &grid_out, double deltat, double time) {

    //grid -> Y1
    //kappa -> weighted_deltat 
    std::array<double, 5> kappa;

    kappa[0] = RKConstants::alpha1 * deltat;
    kappa[1] = RKConstants::alpha2 * deltat;
    kappa[2] = RKConstants::alpha3 * deltat;
    kappa[3] = RKConstants::alpha4 * deltat;
    kappa[4] = RKConstants::alpha5 * deltat;

    
    //BUFFERS
    StaggeredGrid<Addressing_T::STANDARD> Y2(grid.nodes);
    StaggeredGrid<Addressing_T::STANDARD> Y3(grid.nodes);

    #ifdef ForcingTerm 
    ForcingTerm forcingterm(Re);
    #endif

    //Y2.
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {

                #ifdef ForcingTerm
                Y2(Component::U, i, j, k) = grid(Component::U, i, j, k) + rhs(grid, Component::U, i, j, k, kappa[0], Re, forcingterm, time);
                Y2(Component::V, i, j, k) = grid(Component::V, i, j, k) + rhs(grid, Component::V, i, j, k, kappa[0], Re, forcingterm, time);
                Y2(Component::W, i, j, k) = grid(Component::W, i, j, k) + rhs(grid, Component::W, i, j, k, kappa[0], Re, forcingterm, time);
                #else
                Y2(Component::U, i, j, k) = grid(Component::U, i, j, k) + rhs(grid, Component::U, i, j, k, kappa[0], Re);
                Y2(Component::V, i, j, k) = grid(Component::V, i, j, k) + rhs(grid, Component::V, i, j, k, kappa[0], Re);
                Y2(Component::W, i, j, k) = grid(Component::W, i, j, k) + rhs(grid, Component::W, i, j, k, kappa[0], Re);
                #endif

            }
        }
    }

    //Y3.
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {

                #ifdef ForcingTerm
                Y3(Component::U, i, j, k) = Y2(Component::U, i, j, k) + rhs(Y2, Component::U, i, j, k, kappa[1], Re, forcingterm, (time+ Kappa[0])) -
                                            rhs(grid, Component::U, i, j, k, kappa[2], Re, forcingTerm,time);
                Y3(Component::V, i, j, k) = Y2(Component::V, i, j, k) + rhs(Y2, Component::V, i, j, k, kappa[1], Re, forcingterm, (time+ Kappa[0])) -
                                            rhs(grid, Component::V, i, j, k, kappa[2], Re, forcingTerm,time);
                Y3(Component::W, i, j, k) = Y2(Component::W, i, j, k) + rhs(Y2, Component::W, i, j, k, kappa[1], Re, forcingterm, (time+ Kappa[0])) -
                                            rhs(grid, Component::W, i, j, k, kappa[2], Re, forcingTerm,time);
                #else
                Y3(Component::U, i, j, k) = Y2(Component::U, i, j, k) + rhs(Y2, Component::U, i, j, k, kappa[1], Re) -
                                            rhs(grid, Component::U, i, j, k, kappa[2], Re);
                Y3(Component::V, i, j, k) = Y2(Component::V, i, j, k) + rhs(Y2, Component::V, i, j, k, kappa[1], Re) -
                                            rhs(grid, Component::V, i, j, k, kappa[2], Re);
                Y3(Component::W, i, j, k) = Y2(Component::W, i, j, k) + rhs(Y2, Component::W, i, j, k, kappa[1], Re) -
                                            rhs(grid, Component::W, i, j, k, kappa[2], Re);
                #endif



            }
        }
    }

    //u(n+1)    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {

                #ifdef ForcingTerm
                grid_out(Component::U, i, j, k) =
                        Y3(Component::U, i, j, k) + rhs(Y3, Component::U, i, j, k, kappa[3], Re, forcingterm, (time+ Kappa[4])) -
                        rhs(Y2, Component::U, i, j, k, kappa[1], Re, forcingterm, (time+ Kappa[0]));
                grid_out(Component::V, i, j, k) =
                        Y3(Component::V, i, j, k) + rhs(Y3, Component::V, i, j, k, kappa[3], Re, forcingterm, (time+ Kappa[4])) -
                        rhs(Y2, Component::V, i, j, k, kappa[1], Re, forcingterm, (time+ Kappa[0]));
                grid_out(Component::W, i, j, k) =
                        Y3(Component::W, i, j, k) + rhs(Y3, Component::W, i, j, k, kappa[3], Re, forcingterm, (time+ Kappa[4])) -
                        rhs(Y2, Component::W, i, j, k, kappa[1], Re, forcingterm, (time+ Kappa[0]));
                #else
                grid_out(Component::U, i, j, k) =
                        Y3(Component::U, i, j, k) + rhs(Y3, Component::U, i, j, k, kappa[3], Re) -
                        rhs(Y2, Component::U, i, j, k, kappa[1], Re);
                grid_out(Component::V, i, j, k) =
                        Y3(Component::V, i, j, k) + rhs(Y3, Component::V, i, j, k, kappa[3], Re) -
                        rhs(Y2, Component::V, i, j, k, kappa[1], Re);
                grid_out(Component::W, i, j, k) =
                        Y3(Component::W, i, j, k) + rhs(Y3, Component::W, i, j, k, kappa[3], Re) -
                        rhs(Y2, Component::W, i, j, k, kappa[1], Re);
                #endif



            }
        }
    }

}


#endif //AEROHPC_A_RUNGEKUTTA_H
