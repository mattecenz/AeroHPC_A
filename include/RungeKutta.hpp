#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include "data/SolverData.hpp"
#include "equations/MomentumEquation.hpp"
#include "equations/PressureEquation.hpp"

#include "ForcingTerm.hpp"
#include "Boundaries.hpp"

#include "utils/printBuffer.hpp"

// Runge-Kutta method
inline void rungeKutta(const Real time) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const Real t_0 = time + RKConst::beta0 * params.dt;
    const Real t_1 = time + RKConst::beta1 * params.dt;
    const Real t_2 = time + params.dt;


#ifdef ForcingT
    ForcingTerm ft(params.Re, time);
#endif

    Y2star(rkData.buffer_data, rkData.model_data, rkData.model_data)
    apply_boundaries(rkData.buffer_data, t_0, TYPE_VELOCITY);

#ifndef DISABLE_PRESSURE
    P_Eq(params.dt, rkData.buffer_data, rkData.buffer_data)
    apply_boundaries(rkData.buffer_data, t_0, TYPE_PRESSURE);

    Y2(rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.model_data)
    apply_boundaries(rkData.buffer_data, t_0, TYPE_VELOCITY | TYPE_PRESSURE);
#endif

#ifdef ForcingT
    ft.set_time(t_0);
#endif

    Y3star(rkData.model_data, rkData.buffer_data, rkData.buffer_data)
    apply_boundaries(rkData.model_data, t_1, TYPE_VELOCITY);

#ifndef DISABLE_PRESSURE
    P_Eq(params.dt, rkData.model_data, rkData.model_data)
    apply_boundaries(rkData.model_data, t_1, TYPE_PRESSURE);

    Y3(rkData.model_data, rkData.model_data, rkData.model_data, rkData.model_data, rkData.buffer_data)
    apply_boundaries(rkData.model_data, t_1, TYPE_VELOCITY | TYPE_PRESSURE);
#endif


#ifdef ForcingT
    ft.set_time(t_1);
#endif

    U_N1star(rkData.buffer_data, rkData.model_data, rkData.model_data)
    apply_boundaries(rkData.buffer_data, t_2, TYPE_VELOCITY);

#ifndef DISABLE_PRESSURE
    P_Eq(params.dt, rkData.buffer_data, rkData.buffer_data)
    apply_boundaries(rkData.buffer_data, t_2, TYPE_PRESSURE);

    U_N1(rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.model_data);
    apply_boundaries(rkData.buffer_data, t_2, TYPE_VELOCITY | TYPE_PRESSURE);
#endif

    swap(rkData.buffer_data, rkData.model_data);
}

#endif // AEROHPC_A_RUNGEKUTTA_H
