#include <iostream>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include <complex>
#include <cassert>

// Utility function to calculate the squared norm of a 3D vector
double l2_norm_3d(const std::vector<std::vector<std::vector<double>>>& x,
                  const std::vector<std::vector<std::vector<double>>>& x_reconstructed) {
    size_t N = x.size();
    double norm = 0.0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            for (size_t k = 0; k < N; ++k) {
                norm += std::pow(x[i][j][k] - x_reconstructed[i][j][k], 2);
            }
        }
    }
    norm = std::sqrt(norm) / (N * N * N);
    return norm;
}

int main() {
    const int N = 10; // Cube size
    const double L = 1.0; // Length of each dimension
    const double dx = L / N;

    // Step 1: Generate a 3D array x (sine function as an example)
    std::vector<std::vector<std::vector<double>>> x(N,
        std::vector<std::vector<double>>(N, std::vector<double>(N, 0.0)));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                double xi = 2 * M_PI * i / N;
                double yj = 2 * M_PI * j / N;
                double zk = 2 * M_PI * k / N;
                x[i][j][k] = std::sin(xi) * std::sin(yj) * std::sin(zk);
            }
        }
    }

    // Step 2: Initialize b = A * x (here using x directly as b for simplicity)
    std::vector<std::vector<std::vector<double>>> b = x;

    // Step 3: Allocate FFTW complex arrays for 3D
    fftw_complex* b_complex = fftw_alloc_complex(N * N * N);
    fftw_complex* b_tilde = fftw_alloc_complex(N * N * N);

    // Copy data from b into FFTW input array
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                size_t idx = i * N * N + j * N + k;
                b_complex[idx][0] = b[i][j][k]; // Real part
                b_complex[idx][1] = 0.0;       // Imaginary part
            }
        }
    }

    // Step 4: Perform forward FFT
    fftw_plan plan_b = fftw_plan_dft_3d(N, N, N, b_complex, b_tilde, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_b);

    // Step 5: Solve in Fourier space
    std::vector<std::vector<std::vector<double>>> k_squared(N,
        std::vector<std::vector<double>>(N, std::vector<double>(N, 0.0)));
    std::vector<std::complex<double>> x_tilde(N * N * N);

    for (int i = 0; i < N; ++i) {
        int kx = (i < N / 2) ? i : i - N;
        for (int j = 0; j < N; ++j) {
            int ky = (j < N / 2) ? j : j - N;
            for (int k = 0; k < N; ++k) {
                int kz = (k < N / 2) ? k : k - N;

                size_t idx = i * N * N + j * N + k;
                double ksq = - (kx * kx + ky * ky + kz * kz);
                k_squared[i][j][k] = ksq == 0 ? 1e-10 : ksq;

                x_tilde[idx] = std::complex<double>(b_tilde[idx][0], b_tilde[idx][1]) /
                               (-k_squared[i][j][k]);
            }
        }
    }

    // Step 6: Perform inverse FFT
    fftw_complex* x_reconstructed_complex = fftw_alloc_complex(N * N * N);
    for (int i = 0; i < N * N * N; ++i) {
        x_reconstructed_complex[i][0] = x_tilde[i].real();
        x_reconstructed_complex[i][1] = x_tilde[i].imag();
    }

    fftw_plan plan_x = fftw_plan_dft_3d(N, N, N, x_reconstructed_complex, b_complex, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_x);

    // Normalize and store reconstructed x
    std::vector<std::vector<std::vector<double>>> x_reconstructed(N,
        std::vector<std::vector<double>>(N, std::vector<double>(N, 0.0)));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                size_t idx = i * N * N + j * N + k;
                x_reconstructed[i][j][k] = b_complex[idx][0] / (N * N * N); // Normalize
            }
        }
    }

    // Step 7: Compute L2 norm
    double norm = l2_norm_3d(x, x_reconstructed);
    std::cout << "L2 norm: " << norm << std::endl;

    // Print the original and reconstructed values of x pointwise
    std::cout << "Pointwise comparison of original and reconstructed x:" << std::endl;
    for (int i = 5; i < 6; ++i) {
        for (int j = 5; j < 6; ++j) {
            for (int k = 5; k < 6; ++k) {
                std::cout << "x[" << i << "][" << j << "][" << k << "] = "
                        << x[i][j][k] << ", reconstructed = "
                        << x_reconstructed[i][j][k] << std::endl;
            }
        }
    }


    // Cleanup
    fftw_destroy_plan(plan_b);
    fftw_destroy_plan(plan_x);
    fftw_free(b_complex);
    fftw_free(b_tilde);
    fftw_free(x_reconstructed_complex);

    return 0;
}
