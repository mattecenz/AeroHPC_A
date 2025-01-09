#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include "SolverData.hpp"
#include "ForcingTerm.hpp"
#include "Boundaries.hpp"
#include "MomentumEquation.hpp"
#include "PressureEquation.hpp"

#include "printBuffer.hpp"


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
    // TODO APPLY BOUNDARIES ON VELOCITY
    //boundary_cond.apply(model_buff, t_0);


#ifndef DISABLE_PRESSURE
    P_Eq(params.dt, rkData.buffer_data, rkData.buffer_data)
    // TODO APPLY BOUNDARY ON PRESSURE
    // boundary_cond.apply(model_buff, t_0);

    Y2(rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.model_data)
    // TODO APPLY BOUNDARY ON VELOCITY
    // TODO APPLY BOUNDARY ON PRESSURE
    // boundary_cond.apply(model_buff, t_0);
#endif

#ifdef ForcingT
    ft.set_time(t_0);
#endif

    Y3star(rkData.model_data, rkData.buffer_data, rkData.buffer_data)
    //TODO boundary of velocity
    //boundary_cond.apply(model, t_1);

#ifndef DISABLE_PRESSURE
    P_Eq(params.dt, rkData.model_data, rkData.model_data)
    //TODO boundary of pressure
    //boundary_cond.apply(model, t_1);

    Y3(rkData.model_data, rkData.model_data, rkData.model_data, rkData.model_data, rkData.buffer_data)
    //TODO boundary of velocity
    //TODO boundary of pressure
    //boundary_cond.apply(model, t_1);
#endif


#ifdef ForcingT
    ft.set_time(t_1);
#endif

    U_N1star(rkData.buffer_data, rkData.model_data, rkData.model_data)
    //TODO boundary of velocity
    //boundary_cond.apply(model_buff, t_2);

#ifndef DISABLE_PRESSURE
    P_Eq(params.dt, rkData.buffer_data, rkData.buffer_data)
    //TODO boundary of pressure
    //boundary_cond.apply(model_buff, t_2);

    U_N1(rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.model_data);
    //TODO boundary of velocity
    //TODO boundary of pressure
    //boundary_cond.apply(model_buff, t_2);
#endif

    swap(rkData.buffer_data, rkData.model_data);
}

#endif // AEROHPC_A_RUNGEKUTTA_H
