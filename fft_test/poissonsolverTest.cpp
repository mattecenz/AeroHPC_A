#include <iostream>
#include <vector>
#include <cmath>

#include "C2Decomp.hpp"
#include "poissonsolver.hpp"


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

    int N = 120;
    double L = 1.0;
    int pRow = 0, pCol=0;
    bool periodicBC[3] = {true, true, true};
    C2Decomp *c2d = new C2Decomp(N, N, N, pRow, pCol, periodicBC);

    int xsize = c2d->xSize[0];
    int ysize = c2d->xSize[1];
    int zsize = c2d->xSize[2];

    double *b = new double[xsize * ysize * zsize];
    double *X = new double[xsize * ysize * zsize];

    // initialize b as 
    // b = -3 * cos(x) * cos(y) * cos(z)
    for (size_t x = 0; x < xsize; ++x) {
        for (size_t y = 0; y < ysize; ++y) {
            for (size_t z = 0; z < zsize; ++z) {
                size_t index = x * (ysize * zsize) + y * zsize + z; // 3D to 1D mapping
                b[index] = -3.0 * std::cos(M_PI * x) * std::cos(M_PI * y) * std::cos(M_PI * z);
            }
        }
    }

    for (size_t i = 0; i < xsize * ysize * zsize; ++i) {
        X[i] = 0.0;
    }


    poissonSolver solver(N, L, b, c2d);
    solver.solve(X);


    double *X_exact = new double[xsize * ysize * zsize];
    // define X_exact as 
    // X_exact = cos(x) * cos(y) * cos(z)
    for (size_t x = 0; x < xsize; ++x) {
        for (size_t y = 0; y < ysize; ++y) {
            for (size_t z = 0; z < zsize; ++z) {
                size_t index = x * (ysize * zsize) + y * zsize + z; // 3D to 1D mapping
                X_exact[index] = std::cos(M_PI * x) * std::cos(M_PI * y) * std::cos(M_PI * z);
            }
        }
    }

    // Compute the L2 norm
    double local_sum = 0.0;
    for (size_t i = 0; i < xsize * ysize * zsize; ++i) {
        double diff = X[i] - X_exact[i];
        local_sum += diff * diff;
    }

    // Reduce the local sum to compute the global sum
    double global_sum = 0.0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (!mpiRank) {
        double l2_norm = std::sqrt(global_sum / (N * N * N)); // Divide by total number of elements
        cout << "L2 norm between X and X_exact: " << l2_norm << endl;
    }

    if (!mpiRank){
        cout << X[0] << endl;
        cout << X_exact[0] << endl;
    }

    // Cleanup
    delete[] b;
    delete[] X;
    delete[] X_exact;

    // Finalize MPI
    MPI_Finalize();
    return 0;
}