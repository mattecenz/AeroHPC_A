/**
 @file chronoUtils.hpp
 @brief Chronometer Utils
 */

#ifndef AEROHPC_A_CHRONOUTILS_H
#define AEROHPC_A_CHRONOUTILS_H

#include <chrono>

#define timeSpan(start, end) std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
#define chronow() std::chrono::high_resolution_clock::now()

#define measure_t float
#define chrono_start(measure_name) \
                    auto t##measure_name = chronow()
#define chrono_stop(measure_name) \
                    measure_t measure_name = timeSpan(t##measure_name, chronow()) / 1000000000.0

#endif //AEROHPC_A_CHRONOUTILS_H
