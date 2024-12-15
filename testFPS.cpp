#include <iostream>
#include <vector>

#include "C2Decomp.hpp"
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
    int pRow = 0, pCol=0;
    bool periodicBC[3] = {true, true, true};
    C2Decomp *c2d = new C2Decomp(N, N, N, pRow, pCol, periodicBC);

    double *b = new double[N * N * N];
    double *X = new double[N * N * N];

    poissonSolver solver(N, L, b, c2d);
    solver.solve(X);

    // Finalize MPI
    MPI_Finalize();
    return 0;
}