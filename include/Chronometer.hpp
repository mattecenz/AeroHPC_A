#ifndef AEROHPC_A_CHRONOMETER_H
#define AEROHPC_A_CHRONOMETER_H

#include <chrono>
#define timeSpan(start, end) std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
#define chronow() std::chrono::high_resolution_clock::now()

#define measure_t float
#define code_span(code) code
#define chrono_sect(measure_name, code_span)  measure_t measure_name;                        \
                    auto t##measure_name = chronow();                           \
                    code_span                                           \
                    measure_name = timeSpan(t##measure_name, chronow()) / 1000000000.0


#endif //AEROHPC_A_CHRONOMETER_H
