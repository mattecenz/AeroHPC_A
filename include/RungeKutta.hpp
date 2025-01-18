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

    enabledBufferPrinter.print(rkData.buffer_data, "Y2star no bounds", BufferPrinter::PRINT_VELOCITY);
    apply_boundaries(rkData.buffer_data, t_0, TYPE_VELOCITY);
    enabledBufferPrinter.print(rkData.buffer_data, "Y2star", BufferPrinter::PRINT_VELOCITY);

#if !DISABLE_PRESSURE
    P_Eq(consts.inv_k_0, rkData.buffer_data, rkData.buffer_data)

    enabledBufferPrinter.print(rkData.buffer_data, "PHI2-PN incomplete", BufferPrinter::PRINT_PRESSURE);
    apply_boundaries(rkData.buffer_data, t_0, TYPE_PRESSURE);
    enabledBufferPrinter.print(rkData.buffer_data, "PHI2-PN", BufferPrinter::PRINT_PRESSURE);

    Y2(rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.model_data)

    enabledBufferPrinter.print(rkData.buffer_data, "Y2 and PHI 2 no bounds",
                               BufferPrinter::PRINT_VELOCITY | BufferPrinter::PRINT_PRESSURE);
    apply_boundaries(rkData.buffer_data, t_0, TYPE_VELOCITY | TYPE_PRESSURE);
    enabledBufferPrinter.print(rkData.buffer_data, "Y2 and PHI 2 bounds",
                               BufferPrinter::PRINT_VELOCITY | BufferPrinter::PRINT_PRESSURE);
#endif

#if ForcingT
    ft.set_time(t_0);
#endif

    Y3star(rkData.model_data, rkData.buffer_data, rkData.buffer_data)

    enabledBufferPrinter.print(rkData.model_data, "Y3star no bounds", BufferPrinter::PRINT_VELOCITY);
    apply_boundaries(rkData.model_data, t_1, TYPE_VELOCITY);
    enabledBufferPrinter.print(rkData.model_data, "Y3star", BufferPrinter::PRINT_VELOCITY);

#if !DISABLE_PRESSURE
    P_Eq(consts.inv_k_3, rkData.model_data, rkData.model_data)

    enabledBufferPrinter.print(rkData.model_data, "PHI3-PHI2 incomplete", BufferPrinter::PRINT_PRESSURE);
    apply_boundaries(rkData.model_data, t_1, TYPE_PRESSURE);
    enabledBufferPrinter.print(rkData.model_data, "PHI3-PHI2", BufferPrinter::PRINT_PRESSURE);

    Y3(rkData.model_data, rkData.model_data, rkData.model_data, rkData.model_data, rkData.buffer_data)

    enabledBufferPrinter.print(rkData.model_data, "Y3 and PHI 3 no bounds",
                               BufferPrinter::PRINT_VELOCITY | BufferPrinter::PRINT_PRESSURE);
    apply_boundaries(rkData.model_data, t_1, TYPE_VELOCITY | TYPE_PRESSURE);
    enabledBufferPrinter.print(rkData.model_data, "Y3 and PHI 3 bounds",
                               BufferPrinter::PRINT_VELOCITY | BufferPrinter::PRINT_PRESSURE);
#endif


#if ForcingT
    ft.set_time(t_1);
#endif

    U_N1star(rkData.buffer_data, rkData.model_data, rkData.model_data)

    enabledBufferPrinter.print(rkData.buffer_data, "U_N1 no bounds", BufferPrinter::PRINT_VELOCITY);
    apply_boundaries(rkData.buffer_data, t_2, TYPE_VELOCITY);
    enabledBufferPrinter.print(rkData.buffer_data, "U_N1", BufferPrinter::PRINT_VELOCITY);


#if !DISABLE_PRESSURE
    P_Eq(consts.inv_k_6, rkData.buffer_data, rkData.buffer_data)

    enabledBufferPrinter.print(rkData.buffer_data, "PN1-PHI3 incomplete", BufferPrinter::PRINT_PRESSURE);
    apply_boundaries(rkData.buffer_data, t_2, TYPE_PRESSURE);
    enabledBufferPrinter.print(rkData.buffer_data, "PN1-PHI3", BufferPrinter::PRINT_PRESSURE);

    U_N1(rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.buffer_data, rkData.model_data);

    enabledBufferPrinter.print(rkData.buffer_data, "UN1 and PN1 no bounds",
                               BufferPrinter::PRINT_VELOCITY | BufferPrinter::PRINT_PRESSURE);
    apply_boundaries(rkData.buffer_data, t_2, TYPE_VELOCITY | TYPE_PRESSURE);
    enabledBufferPrinter.print(rkData.buffer_data, "UN1 and PN1 bounds",
                               BufferPrinter::PRINT_VELOCITY | BufferPrinter::PRINT_PRESSURE);
#endif

    swap(rkData.buffer_data, rkData.model_data);
}

#endif // AEROHPC_A_RUNGEKUTTA_H
