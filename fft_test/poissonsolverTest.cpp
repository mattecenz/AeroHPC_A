#include <iostream>
#include <vector>
#include <cmath>

#include "C2Decomp.hpp"
#include "poissonSolver.hpp"


using namespace std;

int main(int argc, char *argv[]) {

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
    	cout << " Testing the Fast Poisson Solver " << endl;
    	cout << "-------------------" << endl;
    	cout << endl;
    }

    int N = 2;
    double L = 1.0;
    const double dx = L / N;
    int pRow = 0, pCol=0;
    bool periodicBC[3] = {true, true, true};
    C2Decomp *c2d = new C2Decomp(N, N, N, pRow, pCol, periodicBC);

    // define the global cube
    double *GlobalGrid = new double[N * N * N];
    for (int i=0; i<N; ++i){
        for (int j=0; j<N; ++j){
            for (int k=0; k<N; ++k){
                double xi = i * dx;
                double yj = j * dx;
                double zk = k * dx;
                int index = k * (N * N) + j * N + i;
                GlobalGrid[index] = - 3 * cos(xi) * cos(yj) * cos(zk);
            }
        }
    }

    MPI_Barrier(MPI_Comm MPI_COMM_WORLD);
    // print Global Grid
    if (!mpiRank){
        for (int i=0; i<N; ++i){
            for (int j=0; j<N; ++j){
                for (int k=0; k<N; ++k){
                    int index = k * (N * N) + j * N + i;
                    cout << "i: " << i << ", j: " << j << ", k: " << k << "\t Val = " << GlobalGrid[index] << endl;
                }
            }
        }
    }

    int xsize = c2d->xSize[0];
    int ysize = c2d->xSize[1];
    int zsize = c2d->xSize[2];

    double *b = new double[xsize * ysize * zsize];
    double *X = new double[xsize * ysize * zsize];

    // init b
    int xstart = (mpiRank % (N / xsize)) * xsize;
    int ystart = ((mpiRank / (N / xsize)) % (N / ysize)) * ysize;
    int zstart = (mpiRank / ((N / xsize) * (N / ysize))) * zsize;
    for (int i = 0; i < xsize; ++i) {
        for (int j = 0; j < ysize; ++j) {
            for (int k = 0; k < zsize; ++k) {
                double xi = (xstart + i) * dx;
                double yj = (ystart + j) * dx;
                double zk = (zstart + k) * dx;
                int index = k * (ysize * xsize) + j * xsize + i;
                b[index] = - 3 * cos(xi) * cos(yj) * cos(zk);
            }
        }
    }

    MPI_Barrier(MPI_Comm MPI_COMM_WORLD);
    // print b
    cout << "\n\nprint inital b" << endl;
    for (int rrank = 0; rrank<totRank; rrank++){
        if (mpiRank == rrank){
            for (int i=0; i<xsize; ++i){
                for (int j=0; j<ysize; ++j){
                    for (int k=0; k<zsize; ++k){
                        int index = k * (ysize * xsize) + j * xsize + i;
                        cout << "rank: " << rrank << " i: " << i << ", j: " << j << ", k: " << k << "\t Val = " << b[index] << endl;
                    }
                }
            }
        }
    }


    poissonSolver solver(N, L, c2d);
    solver.setB(b);
    solver.solve(X);
    // solver.performFFT();
    // MPI_Barrier(MPI_Comm MPI_COMM_WORLD);
    // solver.solveEigenvalues();
    // MPI_Barrier(MPI_Comm MPI_COMM_WORLD);
    // solver.performIFFT();
    
    // print b
    cout << "\nprint b after IFFT" << endl;
    MPI_Barrier(MPI_Comm MPI_COMM_WORLD);
    for (int rrank = 0; rrank<totRank; rrank++){
        if (mpiRank == rrank){
            for (int i=0; i<xsize; ++i){
                for (int j=0; j<ysize; ++j){
                    for (int k=0; k<zsize; ++k){
                        int index = k * (ysize * xsize) + j * xsize + i;
                        cout << "rank: " << rrank << " i: " << i << ", j: " << j << ", k: " << k << "\t Val = " << b[index] << endl;
                    }
                }
            }
        }
    }


    // define the exact solution x

    double *exact_x = new double[xsize * ysize * zsize];

    // init b
    for (int i = 0; i < xsize; ++i) {
        for (int j = 0; j < ysize; ++j) {
            for (int k = 0; k < zsize; ++k) {
                double xi = (xstart + i) * dx;
                double yj = (ystart + j) * dx;
                double zk = (zstart + k) * dx;
                int index = k * (ysize * xsize) + j * xsize + i;
                exact_x[index] = cos(xi) * cos(yj) * cos(zk);
            }
        }
    }

    // print exact x
    cout << "\nprint exact x" << endl;
    MPI_Barrier(MPI_Comm MPI_COMM_WORLD);
    for (int rrank = 0; rrank<totRank; rrank++){
        if (mpiRank == rrank){
            for (int i=0; i<xsize; ++i){
                for (int j=0; j<ysize; ++j){
                    for (int k=0; k<zsize; ++k){
                        int index = k * (ysize * xsize) + j * xsize + i;
                        cout << "rank: " << rrank << " i: " << i << ", j: " << j << ", k: " << k << "\t Val = " << exact_x[index] << endl;
                    }
                }
            }
        }
    }



    // solver.solve(X);

    // Cleanup
    delete[] b;
    delete[] X;

    // Finalize MPI
    MPI_Finalize();
    return 0;
}