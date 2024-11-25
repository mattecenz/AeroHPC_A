#include <iostream>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include <complex>
#include <cassert>
#include <algorithm>

using namespace std;

int main() {
    cout << "Poisson Solver using FFT" << endl;

    const int N = 20;  // Number of grid points
    const double L = 1.0;  // Length of the domain
    const double dx = L / N;  // Grid spacing

    // Define the source term b (discretized f(x))
    vector<double> b(N);
    for (int i = 0; i < N; ++i) {
        double x = i * dx;  // Grid point
        b[i] = sin(2 * M_PI * x);  // Example source term: f(x) = sin(2*pi*x)
        // b[0] = 0.0;
    }

    // Allocate memory for FFT input/output
    fftw_complex* b_complex = fftw_alloc_complex(N);
    fftw_complex* b_tilde = fftw_alloc_complex(N);

    // Copy b into b_complex as the real part, with imaginary part initialized to 0
    for (int i = 0; i < N; ++i) {
        b_complex[i][0] = b[i];  // Real part
        b_complex[i][1] = 0.0;   // Imaginary part
    }

    // Step 1: compute b_tilde by performing the Fourier transform of b
    // Perform FFT on b to compute b_tilde
    fftw_plan plan_b = fftw_plan_dft_1d(N, b_complex, b_tilde, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_b);

    // Step 2: compute u_tilde by dividing by the fourier rep (-k^2)
    // Solve for u_tilde in Fourier space
    vector<complex<double>> u_tilde(N);
    vector<double> k(N);  // Wave numbers
    for (int i = 0; i < N; ++i) {
        // Compute wave numbers (k = 2*pi*n/L, with periodicity adjustments)
        k[i] = (i < N / 2) ? (2 * M_PI * i / L) : (2 * M_PI * (i - N) / L);
    }

    for (int i = 0; i < N; ++i) {
        if (i == 0) {
            // Handle k = 0 case (mean component)
            u_tilde[i] = complex<double>(0.0, 0.0);
        } else {
            // Solve u_tilde = b_tilde / (-k^2)
            double k_squared = k[i] * k[i];
            u_tilde[i] = complex<double>(b_tilde[i][0], b_tilde[i][1]) / (-k_squared);
        }
    }

    // Step 3: compute u by performing the inverse fft
    fftw_complex* u_tilde_complex = fftw_alloc_complex(N);
    fftw_complex* u = fftw_alloc_complex(N);

    // Copy u_tilde into u_tilde_complex
    for (int i = 0; i < N; ++i) {
        u_tilde_complex[i][0] = u_tilde[i].real();
        u_tilde_complex[i][1] = u_tilde[i].imag();
    }

    // Inverse FFT to transform back to real space
    fftw_plan plan_u = fftw_plan_dft_1d(N, u_tilde_complex, u, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_u);

    // Normalize the output
    vector<double> u_real(N);
    for (int i = 0; i < N; ++i) {
        u_real[i] = u[i][0] / N;  // Divide by N to normalize
    }

    // Print the solution (u)
    cout << "Solution u(x):" << endl;
    for (int i = 0; i < 10; ++i) {  // Print first 10 values for brevity
        cout << u_real[i] << " ";
    }
    cout << endl;

    // Cleanup
    fftw_destroy_plan(plan_b);
    fftw_destroy_plan(plan_u);
    fftw_free(b_complex);
    fftw_free(b_tilde);
    fftw_free(u_tilde_complex);
    fftw_free(u);

    return 0;
}
