#include <iostream>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include <complex>
#include <cassert>

using namespace std;

int main() {
    const int N = 10;  // Cube size
    const double L = 1.0;  // Length of each dimension
    const double dx = L / N;

    // Step 1: Define the source term b (discretized f(x, y, z))
    vector<vector<vector<double>>> b(N,
        vector<vector<double>>(N, vector<double>(N, 0.0)));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                double xi = i * dx;
                double yj = j * dx;
                double zk = k * dx;
                b[i][j][k] = sin(2 * M_PI * xi) * sin(2 * M_PI * yj) * sin(2 * M_PI * zk);  // Example f(x, y, z)
            }
        }
    }

    // Step 2: Allocate FFTW complex arrays for 3D FFT
    fftw_complex* b_complex = fftw_alloc_complex(N * N * N);
    fftw_complex* b_tilde = fftw_alloc_complex(N * N * N);

    // Copy data from b into FFTW input array
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                size_t idx = i * N * N + j * N + k;
                b_complex[idx][0] = b[i][j][k];  // Real part
                b_complex[idx][1] = 0.0;        // Imaginary part
            }
        }
    }

    // Step 3: Perform forward FFT
    fftw_plan plan_b = fftw_plan_dft_3d(N, N, N, b_complex, b_tilde, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_b);

    // Step 4: Solve in Fourier space
    vector<complex<double>> u_tilde(N * N * N);
    for (int i = 0; i < N; ++i) {
        int kx = (i < N / 2) ? i : i - N;
        for (int j = 0; j < N; ++j) {
            int ky = (j < N / 2) ? j : j - N;
            for (int k = 0; k < N; ++k) {
                int kz = (k < N / 2) ? k : k - N;

                size_t idx = i * N * N + j * N + k;
                double k_squared = kx * kx + ky * ky + kz * kz;
                if (k_squared == 0) k_squared = 1e-10;  // Avoid division by zero

                u_tilde[idx] = complex<double>(b_tilde[idx][0], b_tilde[idx][1]) / (-k_squared);
            }
        }
    }

    // Step 5: Perform inverse FFT
    fftw_complex* u_reconstructed_complex = fftw_alloc_complex(N * N * N);
    for (int i = 0; i < N * N * N; ++i) {
        u_reconstructed_complex[i][0] = u_tilde[i].real();
        u_reconstructed_complex[i][1] = u_tilde[i].imag();
    }

    fftw_plan plan_u = fftw_plan_dft_3d(N, N, N, u_reconstructed_complex, b_complex, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_u);

    // Normalize and store the solution
    vector<vector<vector<double>>> u(N,
        vector<vector<double>>(N, vector<double>(N, 0.0)));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                size_t idx = i * N * N + j * N + k;
                u[i][j][k] = b_complex[idx][0] / (N * N * N);  // Normalize
            }
        }
    }

    // Optional: Output a slice of the solution for inspection
    cout << "Slice of solution u at z=0:" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << u[i][j][0] << " ";
        }
        cout << endl;
    }

    // Cleanup
    fftw_destroy_plan(plan_b);
    fftw_destroy_plan(plan_u);
    fftw_free(b_complex);
    fftw_free(b_tilde);
    fftw_free(u_reconstructed_complex);

    return 0;
}
