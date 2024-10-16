#ifndef AEROHPC_A_TRAITS_H
#define AEROHPC_A_TRAITS_H

/**
 * Enum for easy modification of the class addressing
 */
enum class Addressing_T {
    /**
     * Multidimensional grid flattened into a single dimension array composed by cells of 4 values:
     * [u000, v000, w000, p000, u100, v100, w100, p100,  ... ]
     */
    STANDARD = 0
};

/**
 * Typedef for real values
 */
typedef float Real;

#endif //AEROHPC_A_TRAITS_H
