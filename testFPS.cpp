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
    int pRow = 2, pCol=2;
    bool periodicBC[3] = {true, true, true};
    C2Decomp *c2d = new C2Decomp(N, N, N, pRow, pCol, periodicBC);

    int xsize = c2d->xSize[0];
    int ysize = c2d->xSize[1];
    int zsize = c2d->xSize[2];

    std::cout << "decomp:" << std::endl;
    std::cout << xsize << std::endl;
    std::cout << ysize << std::endl;
    std::cout << zsize << std::endl;

    double *b = new double[xsize * ysize * zsize];
    double *X = new double[xsize * ysize * zsize];

    poissonSolver solver(N, L, b, c2d);
    solver.solve(X);

    // solver.setB(b);

    // Finalize MPI
    MPI_Finalize();
    return 0;
}