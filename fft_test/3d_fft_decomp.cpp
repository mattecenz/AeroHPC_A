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

using namespace std;

#include "C2Decomp.hpp"

 
int main(int argc, char *argv[]){
   int ierr, totRank, mpiRank;

    //Initialize MPI
    ierr = MPI_Init( &argc, &argv);

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


    const int N = 20;  // Cube size
    const double L = 1.0;  // Length of each dimension
    const double dx = L / N;

    // Step 1: Define the source term b (discretized f(x, y, z))
    double *b = (double *)malloc(sizeof(double) * N * N * N);
        
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                double xi = i * dx;
                double yj = j * dx;
                double zk = k * dx;
                // Compute the 1D linear index for the 3D grid
                int index = i * (N * N) + j * N + k;
                // Fill the real and imaginary parts of the FFTW complex array
                // b[index][0] = sin(2 * M_PI * xi) * sin(2 * M_PI * yj) * sin(2 * M_PI * zk);  // Real part
                b[index] = i*N*N + j*N + k;  // Real part
            }
        }
    }

    // if (!mpiRank){
    //     cout << "b:  " << endl;
    //     for (int i=0; i<N*N*N; i++)
    //         cout << b[i] << endl;
    // }

    int pRow = 0, pCol = 0;
    bool periodicBC[3] = {true, true, true};

    if(!mpiRank) cout << "initializing " << endl;
    C2Decomp *c2d;
    c2d = new C2Decomp(N, N, N, pRow, pCol, periodicBC);
    if(!mpiRank) cout << "done initializing " << endl;


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

    // if (!mpiRank){
    //     cout << "u1:  " << endl;
    //     for (int i=0; i<xSize[0]*xSize[1]*xSize[2]; i++)
    //         cout << u1[i] << endl;
    // }

    // Plan for FFT along x-direction
    int local_nx = xSize[0];
    // cout << "local x" << local_nx << endl;
    int batch_size = xSize[1] * xSize[2];

    fftw_plan fft_x_plan = fftw_plan_many_dft(
        1,                      // Rank of the transform (1D FFT)
        &local_nx,              // Size of the transform (local x-size)
        batch_size,             // Number of transforms (y-z grid size)
        reinterpret_cast<fftw_complex *>(u1), // Input array
        NULL,                   // Input strides
        1,                      // Distance between elements in 1 transform
        local_nx,               // Distance between transforms
        reinterpret_cast<fftw_complex *>(u1), // Output array
        NULL,                   // Output strides
        1,                      // Distance between output elements in 1 transform
        local_nx,               // Distance between transforms in output
        FFTW_FORWARD,           // FFT direction
        FFTW_ESTIMATE           // Plan creation mode
    );

    // Execute the FFT along x-direction
    fftw_execute(fft_x_plan);

    // Cleanup the FFTW plan
    fftw_destroy_plan(fft_x_plan);

    // u1=input, u2=output
    c2d->transposeX2Y(u1, u2);
    MPI_Barrier(MPI_COMM_WORLD);

    // if (!mpiRank){
    //     cout << "u2:  " << endl;
    //     for (int i=0; i<ySize[0]*ySize[1]*ySize[2]; i++)
    //         cout << u2[i] << endl;
    // }

    // Plan for FFT along y-direction
    int local_ny = ySize[1];                  // Size of the y-direction (local)
    int batch_size_y = ySize[0] * ySize[2];  // Number of independent FFTs (x-z grid size)

    fftw_plan fft_y_plan = fftw_plan_many_dft(
        1,                          // Rank of the transform (1D FFT)
        &local_ny,                  // Pointer to the size of the transform (y-direction size)
        batch_size_y,               // Number of transforms (x-z grid size)
        reinterpret_cast<fftw_complex *>(u2), // Input array
        NULL,                       // Input strides (default contiguous)
        1,                          // Distance between elements in 1 transform
        local_ny,                   // Distance between transforms
        reinterpret_cast<fftw_complex *>(u2), // Output array (in-place transform)
        NULL,                       // Output strides (default contiguous)
        1,                          // Distance between output elements in 1 transform
        local_ny,                   // Distance between transforms in output
        FFTW_FORWARD,               // FFT direction
        FFTW_ESTIMATE               // Plan creation mode
    );

    // Execute the FFT along y-direction
    fftw_execute(fft_y_plan);

    // Cleanup the FFTW plan
    fftw_destroy_plan(fft_y_plan);

    c2d->transposeY2Z(u2, u3);
    MPI_Barrier(MPI_COMM_WORLD);

    // if (!mpiRank){
    //     cout << "u3:  " << endl;
    //     for (int i=0; i<zSize[0]*zSize[1]*zSize[2]; i++)
    //         cout << u3[i] << endl;
    // }

    // Plan for FFT along z-direction
    int local_nz = zSize[2];                  // Size of the z-direction (local)
    int batch_size_z = zSize[0] * zSize[1];  // Number of independent FFTs (x-y grid size)

    fftw_plan fft_z_plan = fftw_plan_many_dft(
        1,                          // Rank of the transform (1D FFT)
        &local_nz,                  // Pointer to the size of the transform (z-direction size)
        batch_size_z,               // Number of transforms (x-y grid size)
        reinterpret_cast<fftw_complex *>(u3), // Input array
        NULL,                       // Input strides (default contiguous)
        1,                          // Distance between elements in 1 transform
        local_nz,                   // Distance between transforms
        reinterpret_cast<fftw_complex *>(u3), // Output array (in-place transform)
        NULL,                       // Output strides (default contiguous)
        1,                          // Distance between output elements in 1 transform
        local_nz,                   // Distance between transforms in output
        FFTW_FORWARD,               // FFT direction
        FFTW_ESTIMATE               // Plan creation mode
    );

    // Execute the FFT along z-direction
    fftw_execute(fft_z_plan);

    // Cleanup the FFTW plan
    fftw_destroy_plan(fft_z_plan);

    c2d->transposeZ2Y(u3, u2);
    MPI_Barrier(MPI_COMM_WORLD);

    // if (!mpiRank){
    //     cout << "u2:  " << endl;
    //     for (int i=0; i<ySize[0]*ySize[1]*ySize[2]; i++)
    //         cout << u2[i] << endl;
    // }

    c2d->transposeY2X(u2, u1);
    MPI_Barrier(MPI_COMM_WORLD);

    // if (!mpiRank){
    //     cout << "u1:  " << endl;
    //     for (int i=0; i<xSize[0]*xSize[1]*xSize[2]; i++)
    //         cout << u1[i] << endl;
    // }

    // Allocate space for b_new on rank 0
    double *b_new = nullptr;
    if (mpiRank == 0) {
        b_new = (double *)malloc(sizeof(double) * N * N * N);
    }

    // Calculate send counts and displacements for MPI_Gatherv
    std::vector<int> recv_counts(totRank), displacements(totRank);
    int local_size = xSize[0] * xSize[1] * xSize[2]; // Local data size

    // for (int r = 0; r < totRank; ++r) {
    //     // Each rank's contribution is its local xSize product
    //     recv_counts[r] = xSize[0] * xSize[1] * xSize[2];
    //     // if (r == totRank)
    //     // displacements[r] = (c2d->xStart[2] * N * N) + (c2d->xStart[1] * N) + c2d->xStart[0];
    // }

    // displacements[0] = 0;
    // displacements[1] = 4;
    // displacements[2] = 2;
    // displacements[3] = 6;

    // // Gather the decomposed data from all ranks into b_new on rank 0
    // MPI_Gatherv(
    //     u1,                                   // Send buffer (local data on each rank)
    //     local_size,                           // Number of elements to send
    //     MPI_DOUBLE,                           // MPI data type
    //     b_new,                                // Receive buffer (on rank 0)
    //     recv_counts.data(),                   // Number of elements to receive from each rank
    //     displacements.data(),                 // Offsets in the global array for each rank
    //     MPI_DOUBLE,                           // MPI data type
    //     0,                                    // Root rank
    //     MPI_COMM_WORLD                        // MPI communicator
    // );

    // // Optionally process or print b_new on rank 0
    // if (mpiRank == 0) {
    //     // Example: Print some values of b_new for verification
    //     for (int i = 0; i < N * N * N; ++i) {
    //         std::cout << "b_new[" << i << "] = " << b_new[i] << endl;
    //     }
    // }

    // c2d->deallocXYZ(u1);
    // c2d->deallocXYZ(u2);
    // c2d->deallocXYZ(u3);

    free(b);
    if (mpiRank == 0) free(b_new);

    //Now lets kill MPI
    MPI_Finalize();

    return 0;
}