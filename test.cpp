#include <mpi.h>
#include "Traits.hpp"
#include "DataExporter.hpp"
#include "runSolver.cpp"
#include "utils/Logger.hpp"
#include <fstream>

#include "data/SolverData.hpp"

#include "testDomainFunctions.cpp"

int testSolver() {
    const int npy = THIS_WORLD_SIZE;
    const int npz = 1;

    const Real dim_x = 1.0;
    const Real dim_y = 1.0;
    const Real dim_z = 2.0;

    const Real origin_x = 0.0;
    const Real origin_y = 0.0;
    const Real origin_z = -1.0;

    const Real deltaT = 1e-4;
    const Real Re = 5000.0;
    const index_t timeSteps = 1000;

    const bool periodicPressureBC[3] = {false, false, false};

    std::vector<SolverInfo::result_t> results;

    std::vector<index_t> nodes = {
        8, 16, 32, 64
    };

    enabledBufferPrinter.initDir();

    SolverInfo solverInfo{
        false,
        -1,
        "solution.vtk",
        "profile.dat",
        {0.0, 0.0, 0.0}
    };

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

        runSolver(solverInfo);

        destroySolverData();

        enabledLogger.printTitle("Data destroyed")
                .closeSection().empty();

        results.push_back(solverInfo.results);
    }

    if (IS_MAIN_PROC) {
        std::ofstream csvFile("output.csv");

        // Print header of csv
        auto header_iter = results.front().begin();
        csvFile << header_iter->first;
        ++header_iter;
        for (; header_iter != results.front().end(); ++header_iter) {
            csvFile << "," << header_iter->first;
        }
        csvFile << std::endl;

        for (const auto &result: results) {
            auto res_iter = result.begin();
            csvFile << res_iter->second;
            ++res_iter;
            for (; res_iter != result.end(); ++res_iter) {
                csvFile << "," << res_iter->second;
            }
            csvFile << std::endl;
        }
    }

    return 0;
}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &THIS_PROC_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &THIS_WORLD_SIZE);

    const int ris = testSolver();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return ris;
}
