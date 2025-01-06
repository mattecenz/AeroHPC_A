#include "Traits.hpp"

Vector operator*(Real b, const Vector &a) {
    return Vector{
            a[0] * b,
            a[1] * b,
            a[2] * b
    };
}

Vector operator*(const Vector &a, Real b) {
    return Vector{
            a[0] * b,
            a[1] * b,
            a[2] * b
    };
}

Real operator*(const Vector &a, const Vector &b) {
    return a[0] * b[0]
           + a[1] * b[1]
           + a[2] * b[2];
}

Vector operator+(const Vector &a, const Vector &b) {
    return Vector{
            a[0] + b[0],
            a[1] + b[1],
            a[2] + b[2]
    };
}

Vector operator-(const Vector &a) {
    return Vector{
            -a[0],
            -a[1],
            -a[2]
    };
}

Vector operator-(const Vector &a, const Vector &b) {
    return Vector{
            a[0] - b[0],
            a[1] - b[1],
            a[2] - b[2]
    };
}