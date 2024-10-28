#ifndef AEROHPC_A_CHRONOMETER_H
#define AEROHPC_A_CHRONOMETER_H

#include <chrono>
#define timeSpan(start, end) std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
#define chronow() std::chrono::high_resolution_clock::now()

#define measure_t float
#define code_span(code) code
#define measure(name, code_span)  measure_t name;                        \
                    auto t##name = chronow();                           \
                    code_span                                           \
                    name = timeSpan(t##name,chronow()) / 1000000000.0


#endif //AEROHPC_A_CHRONOMETER_H
