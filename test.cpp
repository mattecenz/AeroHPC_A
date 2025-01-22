#include <mpi.h>
#include "Traits.hpp"
#include "DataExporter.hpp"
#include "runSolver.cpp"
#include "utils/Logger.hpp"
#include <fstream>

#include "data/SolverData.hpp"

#include "testDomainFunctions.cpp"

int testSolver() {
    const int npy = 1;
    const int npz = 1;

    const Real dim_x = 1.0;
    const Real dim_y = 1.0;
    const Real dim_z = 1.0;

    const Real origin_x = 0.0;
    const Real origin_y = 0.0;
    const Real origin_z = 0.0;

    const Real deltaT = 1e-4;
    const Real Re = 1000.0;
    const index_t timeSteps = 100;

    const bool periodicPressureBC[3] = {false, false, false};

    std::vector<Real> error;


    std::vector<index_t> nodes = {
        8
    };

    enabledBufferPrinter.initDir();

    // wrt dim
    for (index_t n: nodes) {
        enabledLogger.openSection("TEST");

        initSolverData(dim_x, dim_y, dim_z,
                 origin_x, origin_y, origin_z,
                 deltaT, timeSteps, Re,
                 n, n, n,
                 npy, npz,
                 periodicPressureBC,
                 &testDomainData);

        enabledLogger.printTitle("Data initialized");

        Real l2norm = runSolver(0.0, 0.0, 0.0);

        destroySolverData();

        enabledLogger.printTitle("Data destroyed")
                .closeSection().empty();

        error.push_back(l2norm);
    }

    if (IS_MAIN_PROC) {
        std::ofstream csvFile("output.csv");
        csvFile << "step,error" << std::endl;
        for (int i = 0; i < nodes.size(); ++i) csvFile << nodes[i] << "," << error[i] << std::endl;
    }

    return 0;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &THIS_PROC_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &THIS_WORLD_SIZE);

    int ris = testSolver();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return ris;
}
