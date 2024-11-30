#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <cassert>
#include <algorithm>
#include <fstream>

// Function to compute the L2 norm between two arrays
double l2_norm(const double* x, const double* x_reconstructed, int N) {
    double norm = 0.0;
    for (int i = 0; i < N; ++i) {
        norm += std::pow(x[i] - x_reconstructed[i], 2);
    }
    norm = std::sqrt(norm);
    return norm;
}

int main() {
    const int N = 100;  // Size of the array
    const double L = 2.0;  // Length of the domain

    // Step 1: Generate the array x (sine function)
    double* x = fftw_alloc_real(N);
    for (int i = 0; i < N; ++i) {
        double value = (M_PI * i) / (N - 1);  // Equally spaced array from 0 to pi
        x[i] = std::cos(value * 3);
    }

    std::ofstream file_original("x_original.csv");
    if (file_original.is_open()) {
        for (int i = 0; i < N; ++i) {
            file_original << i << "," << x[i] << "\n"; // Write index and value
        }
        file_original.close();
        std::cout << "Original data exported to x_original.csv\n";
    } else {
        std::cerr << "Error: Unable to open file for writing.\n";
    }

    // Step 2: Construct the matrix A
    double** A = new double*[N];
    for (int i = 0; i < N; ++i) {
        A[i] = new double[N]();
    }
    for (int i = 0; i < N; ++i) {
        A[i][i] = -2.0;
    }
    for (int i = 0; i < N - 1; ++i) {
        A[i][i + 1] = 1.0;  // Upper diagonal
        A[i + 1][i] = 1.0;  // Lower diagonal
    }
    A[0][0] = -1.0;           // Neumann boundary condition
    A[N - 1][N - 1] = -1.0;

    // Step 3: Compute b = A * x
    double* b = fftw_alloc_real(N);
    std::fill(b, b + N, 0.0); // Initialize with zeros
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            b[i] += A[i][j] * x[j];
        }
    }

    // Step 4: Apply FFT to b
    double* b_tilde = fftw_alloc_real(N);
    double* result = fftw_alloc_real(N);
    fftw_plan plan_b = fftw_plan_r2r_1d(N, b, b_tilde, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(plan_b);

    // Step 5: Solve for x_tilde in Fourier space
    double* x_tilde = fftw_alloc_real(N);
    double* k = fftw_alloc_real(N);
    for (int i = 0; i < N; ++i) {
        k[i] = (i < N / 2) ? (std::sqrt(2) * M_PI * i / N) : (std::sqrt(2) * M_PI * (i - N) / N);
    }
    k[0] = 1e-10;  // Avoid division by zero
    for (int i = 0; i < N; ++i) {
        if (i != 0) {
            x_tilde[i] = b_tilde[i] / (-std::pow(k[i], 2));
        } else {
            x_tilde[i] = 0.0;
        }
    }

    // Step 6: Perform inverse FFT to get reconstructed x
    double* x_reconstructed = fftw_alloc_real(N);
    fftw_plan plan_x = fftw_plan_r2r_1d(N, x_tilde, x_reconstructed, FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_execute(plan_x);

    for (int i=0; i<N; i++){
        x_reconstructed[i] /= N;
    }

    // Step 7: Compute the L2 norm between x and x_reconstructed
    double norm = l2_norm(x, x_reconstructed, N);
    std::cout << "L2 norm: " << norm << std::endl;

    // Export x_reconstructed to a CSV file



    std::ofstream file("x_reconstructed.csv");
    if (file.is_open()) {
        for (int i = 0; i < N; ++i) {
            file << i << "," << x_reconstructed[i] << "\n"; // Write index and value
        }
        file.close();
        std::cout << "Reconstructed data exported to x_reconstructed.csv\n";
    } else {
        std::cerr << "Error: Unable to open file for writing.\n";
    }

    // Cleanup
    fftw_destroy_plan(plan_b);
    fftw_destroy_plan(plan_x);
    fftw_free(x);
    fftw_free(b);
    fftw_free(b_tilde);
    fftw_free(x_tilde);
    fftw_free(x_reconstructed);
    fftw_free(k);
    for (int i = 0; i < N; ++i) {
        delete[] A[i];
    }
    delete[] A;

    return 0;
}
