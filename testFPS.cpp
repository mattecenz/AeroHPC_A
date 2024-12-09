#include <iostream>
#include <vector>

#include "PoissonSolver.hpp"

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

    int N = 200;
    double L = 1.0;
    double *b = new double[N * N * N];
    double *X = new double[N * N * N];
    bool periodicBC[3] = {true, true, true};

    poissonSolver solver(N, L, b, periodicBC);
    // solver.solve(X);


    return 0;
}