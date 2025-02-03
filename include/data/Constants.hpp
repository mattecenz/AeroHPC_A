#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

struct RKConst {
    static constexpr Real alpha0 = 64.0 / 120.0;
    static constexpr Real alpha1 = 34.0 / 120.0;
    static constexpr Real alpha2 = 50.0 / 120.0;
    static constexpr Real alpha3 = alpha2 - alpha1;
    static constexpr Real alpha4 = 50.0 / 120.0;
    static constexpr Real alpha5 = 90.0 / 120.0;
    static constexpr Real alpha6 = alpha5 - alpha4;
    static constexpr Real beta0 = 64.0 / 120.0;
    static constexpr Real beta1 = 80.0 / 120.0;
};

class Constants {
public:
    Real nu,
            k_0, k_1, k_2, k_3, k_4, k_5, k_6,
            inv_k_0, inv_k_3, inv_k_6;

    Constants(const Real Re, const Real dt) {
        nu = (real(1.0) / Re);

        k_0 = RKConst::alpha0 * dt;
        k_1 = RKConst::alpha1 * dt;
        k_2 = RKConst::alpha2 * dt;
        k_3 = RKConst::alpha3 * dt;
        k_4 = RKConst::alpha4 * dt;
        k_5 = RKConst::alpha5 * dt;
        k_6 = RKConst::alpha6 * dt;
        inv_k_0 = real(1.0) / k_0;
        inv_k_3 = real(1.0) / k_3;
        inv_k_6 = real(1.0) / k_6;
    }
};

#endif //CONSTANTS_HPP
