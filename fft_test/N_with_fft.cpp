#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <fftw3.h>
#include <complex>
#include <cassert>
#include <chrono>

using namespace std;

#include "C2Decomp.hpp"

 
int main(int argc, char *argv[]){

    int ierr, totRank, mpiRank;

    //Initialize MPI
    ierr = MPI_Init(&argc, &argv);

    //Get the number of processes
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &totRank);

    //Get the local rank
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if(!mpiRank){
        cout << endl;
        cout << "-------------------" << endl;
    	cout << " Parallel FFT " << endl;
    	cout << "-------------------" << endl;
    	cout << endl;
    }

    const int N = 100;  // Cube size
    const double L = 1.0;  // Length of each dimension
    const double dx = L / N;

    double *b = (double *)malloc(sizeof(double) * N * N * N);
        
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                double xi = i * dx;
                double yj = j * dx;
                double zk = k * dx;
                // Compute the 1D linear index for the 3D grid
                int index = i * (N * N) + j * N + k;
                b[index] = sin(2 * M_PI * xi) * sin(2 * M_PI * yj) * sin(2 * M_PI * zk);

            }
        }
    }

    int pRow = 0, pCol = 0;
    bool periodicBC[3] = {true, true, true};

    auto start_pre = std::chrono::high_resolution_clock::now();
    if(!mpiRank) cout << "initializing " << endl;
    C2Decomp *c2d;
    c2d = new C2Decomp(N, N, N, pRow, pCol, periodicBC);
    if(!mpiRank) cout << "done initializing " << endl;
    auto end_pre = std::chrono::high_resolution_clock::now();

    if(!mpiRank) cout << "Grid dim: " << N << endl;

    auto start = std::chrono::high_resolution_clock::now();

    double xSize[3], ySize[3], zSize[3];
    xSize[0] = c2d->xSize[0];
    xSize[1] = c2d->xSize[1];
    xSize[2] = c2d->xSize[2];
    ySize[0] = c2d->ySize[0];
    ySize[1] = c2d->ySize[1];
    ySize[2] = c2d->ySize[2];
    zSize[0] = c2d->zSize[0];
    zSize[1] = c2d->zSize[1];
    zSize[2] = c2d->zSize[2];

    double *u1, *u2, *u3;

    c2d->allocX(u1);
    c2d->allocY(u2);
    c2d->allocZ(u3);

    // Distribute the data across the ranks
    for(int kp = 0; kp < xSize[2]; kp++){
        for(int jp = 0; jp < xSize[1]; jp++){
            for(int ip = 0; ip < xSize[0]; ip++){
                int ii = kp*xSize[1]*xSize[0] + jp*xSize[0] + ip;
                
                int global_index = (c2d->xStart[2] + kp) * (N * N) 
                 + (c2d->xStart[1] + jp) * N 
                 + (c2d->xStart[0] + ip);
        
                u1[ii] = b[global_index];  // Assuming you need only the real part here
            }
        }
    }

    // Perform FFT along the X direction (on u1)
    fftw_complex *fft_input, *fft_output;
    fft_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * xSize[0] * xSize[1] * xSize[2]);
    fft_output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * xSize[0] * xSize[1] * xSize[2]);

    // Plan for FFT along x-direction
    fftw_plan plan_x = fftw_plan_r2r_1d(xSize[0], fft_input, fft_output, FFTW_REDFT10, FFTW_ESTIMATE);

    // Perform FFT along x (for each y and z slice)
    for (int jp = 0; jp < xSize[1]; ++jp) {
        for (int kp = 0; kp < xSize[2]; ++kp) {
            for (int ip = 0; ip < xSize[0]; ++ip) {
                int ii = kp * xSize[1] * xSize[0] + jp * xSize[0] + ip;
                fft_input[ip][0] = u1[ii]; // Real part
                fft_input[ip][1] = 0;     // Imaginary part (initially 0)
            }
            fftw_execute(plan_x);  // Execute FFT on each slice

            // Store the FFT result back into u1 (you may store real and imaginary parts here)
            for (int ip = 0; ip < xSize[0]; ++ip) {
                int ii = kp * xSize[1] * xSize[0] + jp * xSize[0] + ip;
                u1[ii] = fft_output[ip][0]; // Storing the real part
            }
        }
    }
    // Transpose u1 to u2 (y-direction)
    c2d->transposeX2Y(u1, u2);
    MPI_Barrier(MPI_COMM_WORLD);

    // Perform FFT along the Y direction (on u2)
    fftw_complex *fft_input_y, *fft_output_y;
    fft_input_y = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ySize[0] * ySize[1] * ySize[2]);
    fft_output_y = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * ySize[0] * ySize[1] * ySize[2]);

    // Plan for FFT along y-direction
    fftw_plan plan_y = fftw_plan_r2r_1d(ySize[0], fft_input_y, fft_output_y, FFTW_REDFT10, FFTW_ESTIMATE);

    // Perform FFT along y (for each x and z slice)
    for (int kp = 0; kp < ySize[2]; ++kp) {
        for (int ip = 0; ip < ySize[0]; ++ip) {
            for (int jp = 0; jp < ySize[1]; ++jp) {
                int ii = kp * ySize[1] * ySize[0] + jp * ySize[0] + ip;
                fft_input_y[ip][0] = u2[ii]; // Real part
                fft_input_y[ip][1] = 0;     // Imaginary part
            }
            fftw_execute(plan_y);  // Execute FFT on each slice

            // Store the FFT result back into u2
            for (int jp = 0; jp < ySize[1]; ++jp) {
                int ii = kp * ySize[1] * ySize[0] + jp * ySize[0] + ip;
                u2[ii] = fft_output_y[ip][0]; // Storing real part
            }
        }
    }

    // Transpose u2 to u3 (z-direction)
    c2d->transposeY2Z(u2, u3);
    MPI_Barrier(MPI_COMM_WORLD);

    // Perform FFT along the Z direction (on u3)
    fftw_complex *fft_input_z, *fft_output_z;
    fft_input_z = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * zSize[0] * zSize[1] * zSize[2]);
    fft_output_z = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * zSize[0] * zSize[1] * zSize[2]);

    // Plan for FFT along z-direction
    fftw_plan plan_z = fftw_plan_r2r_1d(zSize[0], fft_input_z, fft_output_z, FFTW_REDFT10, FFTW_ESTIMATE);

    // Perform FFT along z (for each x and y slice)
    for (int jp = 0; jp < zSize[1]; ++jp) {
        for (int ip = 0; ip < zSize[0]; ++ip) {
            for (int kp = 0; kp < zSize[2]; ++kp) {
                int ii = kp * zSize[1] * zSize[0] + jp * zSize[0] + ip;
                fft_input_z[ip][0] = u3[ii]; // Real part
                fft_input_z[ip][1] = 0;     // Imaginary part
            }
            fftw_execute(plan_z);  // Execute FFT on each slice

            // Store the FFT result back into u3
            for (int kp = 0; kp < zSize[0]; ++kp) {
                int ii = kp * zSize[1] * zSize[0] + jp * zSize[0] + ip;
                u3[ii] = fft_output_z[ip][0]; // Storing the real part
            }
        }
    }

    // Transpose u3 to u2 (reverse z -> y)
    c2d->transposeZ2Y(u3, u2);
    c2d->transposeY2X(u2, u1);


    fftw_free(fft_input);
    fftw_free(fft_output);
    fftw_free(fft_input_y);
    fftw_free(fft_output_y);
    fftw_free(fft_input_z);
    fftw_free(fft_output_z);

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
