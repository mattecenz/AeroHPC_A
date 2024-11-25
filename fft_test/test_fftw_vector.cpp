#include <iostream>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include <complex>
#include <cassert>
#include <algorithm>

// Function to compute the L2 norm between two vectors
double l2_norm(const std::vector<double>& x, const std::vector<double>& x_reconstructed) {
    double norm = 0.0;
    assert(x.size() == x_reconstructed.size());
    for (size_t i = 0; i < x.size(); ++i) {
        norm += std::pow(x[i] - x_reconstructed[i], 2);
    }
    
    norm = std::sqrt(norm) / N;
    return norm;
}

int main() {
    const int N = 1000;  // Size of the array
    const double L = 1.0;  // Length of the domain

    // Step 1: Generate the array x (sine function)
    std::vector<double> x(N);
    for (int i = 0; i < N; ++i) {
        // Create an equally spaced array from 0 to 2 * pi
        double value = (2 * M_PI * i) / (N - 1);  // Note: i goes from 0 to N-1
        x[i] = std::sin(value);
    }

    // print x
    // std::cout << "Original x: ";
    // for (const auto& val : x) {
    //     std::cout << val << ", ";
    // }

    // Step 2: Construct the matrix A with zeros and set the main diagonal to -2
    std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0)); // Initialize with 0s
    for (int i = 0; i < N; ++i) {
        A[i][i] = -2.0;  // Set the main diagonal to -2
    }

    // Set the upper and lower diagonals
    for (int i = 0; i < N - 1; ++i) {
        A[i][i + 1] = 1.0;  // Upper diagonal
        A[i + 1][i] = 1.0;  // Lower diagonal
    }

    // Implement periodicity
    A[0][N - 1] = 1.0;  // (upper-right corner)
    A[N - 1][0] = 1.0;  // (lower-left corner)


    // Printing the matrix A
    // std::cout << "matrix A" << std::endl;
    // for (int i = 0; i < N; ++i) {
    //     for (int j = 0; j < N; ++j) {
    //         std::cout << A[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }


    // Step 3: Compute b = A * x
    std::vector<double> b(N, 0.0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            b[i] += A[i][j] * x[j];
        }
    }

    // std::cout << "b: ";
    // for (const auto& val : b) {
    //     std::cout << val << " ";
    // }

    // Step 3.2: Copy the content of b into a complex array where the imaginary part is initialized to 0
    fftw_complex* b_complex = fftw_alloc_complex(N);
    for (int i = 0; i < N; ++i) {
        b_complex[i][0] = b[i];  // Real part from b
        b_complex[i][1] = 0.0;   // Imaginary part is initialized to 0
    }

    // Step 4: Apply FFT to b to get b_tilde
    fftw_complex* b_tilde = fftw_alloc_complex(N);
    fftw_plan plan_b = fftw_plan_dft_1d(N, b_complex, b_tilde, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_b);

    // Print the values in b_tilde (complex array)
    // std::cout << std::endl;
    // for (int i = 0; i < N; ++i) {   
    // std::cout << "b_tilde[" << i << "] = (" 
    //           << b_tilde[i][0] << ", "   // Real part
    //           << b_tilde[i][1] << ")"   // Imaginary part
    //           << std::endl;
    // }


    // Step 5: Solve for x_tilde in Fourier space by dividing b_tilde by -k^2
    std::vector<std::complex<double>> x_tilde(N);
    std::vector<double> k(N);
    for (int i = 0; i < N; ++i) {
        k[i] = (i < N / 2) ? (2 * M_PI * i / N) : (2 * M_PI * (i - N) / N);  // Use M_PI from cmath
    }
    k[0] = 1e-10;  // Avoid division by zero for the zero frequency component

    // std::cout << std::endl;
    // std::cout << "k: ";
    // for (const auto& val : k) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;

    for (int i = 0; i < N; ++i) {
        if (i != 0) {  // Skip the zero frequency component
            x_tilde[i] = std::complex<double>(b_tilde[i][0], b_tilde[i][1]) / (-std::pow(k[i], 2));
        } else {
            x_tilde[i] = std::complex<double>(0.0, 0.0);  // Set zero frequency component to 0
        }
    }

    // std::cout << "x_tilde: ";
    // for (const auto& val : x_tilde) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;

    // Step 6: Perform inverse FFT on x_tilde to get the reconstructed x
    std::vector<double> x_reconstructed(N, 0.0);
    fftw_complex* x_tilde_complex = fftw_alloc_complex(N);
    for (int i = 0; i < N; ++i) {
        x_tilde_complex[i][0] = x_tilde[i].real();
        x_tilde_complex[i][1] = x_tilde[i].imag();
    }

    fftw_plan plan_x = fftw_plan_dft_c2r_1d(N, x_tilde_complex, x_reconstructed.data(), FFTW_ESTIMATE);
    fftw_execute(plan_x);

    // Normalize the result (since FFTW works with complex numbers)
    for (int i = 0; i < N; ++i) {
        x_reconstructed[i] /= N;
    }

    // Step 7: Compare the original x and the reconstructed x
    double norm = l2_norm(x, x_reconstructed);
    std::cout << "L2 norm: " << norm << std::endl;

    // // Optionally, you can print the arrays and compare
    // std::cout << "Original x: ";
    // for (const auto& val : x) {
    //     std::cout << val << " ";
    // }
    // std::cout << "\nReconstructed x: ";
    // for (const auto& val : x_reconstructed) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;

    // Cleanup
    fftw_destroy_plan(plan_b);
    fftw_destroy_plan(plan_x);
    fftw_free(b_tilde);
    fftw_free(x_tilde_complex);

    return 0;
}
