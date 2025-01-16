#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include "data/SolverData.hpp"
#include "equations/MomentumEquation.hpp"
#include "equations/PressureEquation.hpp"

#include "ForcingTerm.hpp"
#include "Boundaries.hpp"

#include "utils/printBuffer.hpp"

inline void rungeKutta(const Real time) {
    const Real t_0 = time + RKConst::beta0 * params.dt;
    const Real t_1 = time + RKConst::beta1 * params.dt;
    const Real t_2 = time + params.dt;

#if ForcingT
    ForcingTerm ft(params.Re, time);
#endif

    Y2star(rkData.buffer_data, rkData.model_data, rkData.model_data)
    apply_boundaries(rkData.buffer_data, t_0, TYPE_VELOCITY);

#if !DISABLE_PRESSURE
    //CHANGEING 1e-6 to the real coefficient (consts.inv_k_0) makes everything explode
    P_Eq(consts.inv_k_0, rkData.buffer_data, rkData.buffer_data)
    // apply_boundaries(rkData.buffer_data, t_0, TYPE_PRESSURE);

    Y2(rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.model_data)
    apply_boundaries(rkData.buffer_data, t_0, TYPE_VELOCITY /*| TYPE_PRESSURE*/);
#endif

#if ForcingT
    ft.set_time(t_0);
#endif

    Y3star(rkData.model_data, rkData.buffer_data, rkData.buffer_data)
    apply_boundaries(rkData.model_data, t_1, TYPE_VELOCITY);

#if !DISABLE_PRESSURE
    //CHANGEING 1e-6 to the real coefficient (consts.inv_k_3) makes everything explode
    P_Eq(consts.inv_k_3, rkData.model_data, rkData.model_data)
    // apply_boundaries(rkData.model_data, t_1, TYPE_PRESSURE);

    Y3(rkData.model_data, rkData.model_data, rkData.model_data, rkData.model_data, rkData.buffer_data)
    apply_boundaries(rkData.model_data, t_1, TYPE_VELOCITY /*| TYPE_PRESSURE*/);
#endif


#if ForcingT
    ft.set_time(t_1);
#endif

    U_N1star(rkData.buffer_data, rkData.model_data, rkData.model_data)
    apply_boundaries(rkData.buffer_data, t_2, TYPE_VELOCITY);

#if !DISABLE_PRESSURE
    //CHANGEING 1e-6 to the real coefficient (consts.inv_k_6) makes everything explode
    P_Eq(consts.inv_k_6, rkData.buffer_data, rkData.buffer_data)
    // apply_boundaries(rkData.buffer_data, t_2, TYPE_PRESSURE);

    U_N1(rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.model_data);
    apply_boundaries(rkData.buffer_data, t_2, TYPE_VELOCITY /*| TYPE_PRESSURE*/);
#endif

    swap(rkData.buffer_data, rkData.model_data);
}

#endif // AEROHPC_A_RUNGEKUTTA_H
