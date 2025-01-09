#include <mpi.h>
#include "Traits.hpp"
#include "DataExporter.hpp"
#include "runSolver.cpp"
#include <fstream>

#include "SolverData.hpp"

int testSolver(const int rank, const int size) {

    Real c =

    const int npy = 1;
    const int npz = 1;

    const Real dim_x = 1.0;
    const Real dim_y = 1.0;
    const Real dim_z = 1.0;

    const Real origin_x = 0.0;
    const Real origin_y = 0.0;
    const Real origin_z = 0.0;

    const Real deltaT = 10e-4;
    const Real Re = 1000;
    const index_t timeSteps = 1000;

    const bool periodicPressureBC[3] = {true, true, true};

    std::vector<Real> error;

    const boundaryDomainFunctions boundaryDomainFunctions = {
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w},
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w},
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w},
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w},
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w},
        boundaryFaceFunctions{ExactSolution::u, ExactSolution::v, ExactSolution::w}
    };


    std::vector<index_t> nodes = {
        4, 8, 16
    };


    // wrt dim
    for (index_t n : nodes) {
        initData(dim_x, dim_y, dim_z,
                 origin_x, origin_y, origin_z,
                 deltaT, timeSteps, Re,
                 n, n, n,
                 npy, npz,
                 periodicPressureBC);

        Real l2norm = runSolver(rank, size,
                                boundaryDomainFunctions,
                                0.0, 0.0, 0.0);

        destroyData();

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
