#include <mpi.h>
#include "Traits.hpp"
#include "DataExporter.hpp"
#include "runSolver.cpp"
#include <fstream>

int testSolver(const int rank, const int size) {
    const int npy = 1;
    const int npz = 1;

    // dividing the timestep size to half
    const std::vector<Real> deltaTs = {0.001, 0.0005, 0.00025};
    const std::vector<int> nodes = {
        4,8,16 ,32,64
    };
    const Real dim_x = 1.0;
    const Real dim_y = 1.0;
    const Real dim_z = 1.0;
    const Real Re = 1000;
    const index_t timeSteps = 1000;


    std::vector<Real> error;

    // wrt deltaT
    // for (size_t i=0; i<deltaTs.size(); ++i){
    //     Real deltaT = deltaTs[i];
    //     index_t dim = dims[3];
    //     error.push_back(testSolver(deltaT, dim));
    // }

    const boundaryDomainFunctions boundaryDomainFunctions = {
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w},
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w},
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w},
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w},
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w},
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w}
    };

    // wrt dim
    for (int n: nodes) {
        Real deltaT = deltaTs[0]; // first
        Real l2norm = runSolver(rank, size,
                                npy, npz,
                                n, n, n,
                                dim_x, dim_y, dim_z,
                                deltaT, timeSteps, Re,
                                boundaryDomainFunctions,
                                0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0);
        error.push_back(l2norm);
    }


    // // wrt both
    // for (size_t i=0; i<dims.size(); ++i){
    //     Real deltaT = deltaT[i];
    //     index_t dim = dims[i];
    //     error.push_back(testSolver(deltaT, dim));
    // }

    if (!rank) {
        std::ofstream csvFile("output.csv");
        csvFile << "step,error" << std::endl;
        for (int i = 0; i < nodes.size(); ++i) csvFile << nodes[i] << "," << error[i] << std::endl;
    }

    return 0;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int ris = testSolver(rank, size);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return ris;
}
