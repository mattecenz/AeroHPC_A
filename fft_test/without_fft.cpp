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

    const int N = 2;  // Cube size
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
                b[index] = i*N*N + j*N + k; 
            }
        }
    }

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

    // Transpose steps
    c2d->transposeX2Y(u1, u2);
    MPI_Barrier(MPI_COMM_WORLD);
    
    c2d->transposeY2Z(u2, u3);
    MPI_Barrier(MPI_COMM_WORLD);

    c2d->transposeZ2Y(u3, u2);
    MPI_Barrier(MPI_COMM_WORLD);

    c2d->transposeY2X(u2, u1);
    MPI_Barrier(MPI_COMM_WORLD);

    // Now gather all u1 arrays from all processes into one array on the root process (rank 0)
    double *b_new = NULL;

    if(mpiRank == 0) {
        // Allocate space for the full b_new array
        b_new = (double *)malloc(sizeof(double) * N * N * N);
    }

    // Gather the u1 data from each process into the root process
    MPI_Gather(u1, xSize[0] * xSize[1] * xSize[2], MPI_DOUBLE, 
               b_new, xSize[0] * xSize[1] * xSize[2], MPI_DOUBLE, 
               0, MPI_COMM_WORLD);

    // Optionally, print the gathered data on the root process
    if(mpiRank == 0) {
        cout << "b_new (gathered data):" << endl;
        for(int i = 0; i < N * N * N; ++i) {
            cout << b_new[i] << endl;
        }

        // Don't forget to free the memory allocated for b_new
        free(b_new);
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
