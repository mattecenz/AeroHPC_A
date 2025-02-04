#include <mpi.h>
#include "Traits.hpp"
#include "DataExporter.hpp"
#include "runSolver.cpp"
#include "testcaseDomainFunctions.cpp"
#include "utils/Logger.hpp"
#include "data/SolverData.hpp"
#include "data/DomainData.hpp"
#include <fstream>

using namespace std;

int runTestCase(int argc, char **argv) {
    int npy, npz;
    int nx, ny, nz;
    int timeSteps;
    Real deltaT;
    int testCase;

    // Processor 0 read the file
    if (IS_MAIN_PROC) {
        if (argc < 2) {
            cout << "Test case file not found: Abort" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        ifstream inputFile(argv[1]);

        // Check if the file is open
        if (!inputFile.is_open()) {
            std::cerr << "Error: Unable to open the file." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        inputFile >> timeSteps;
        inputFile >> deltaT;
        inputFile >> nx >> ny >> nz;
        inputFile >> npy >> npz;
        inputFile >> testCase;

        inputFile.close();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Then distribute results
    MPI_Bcast(&npy, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&npz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&timeSteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&deltaT, 1, Real_MPI, 0, MPI_COMM_WORLD);
    MPI_Bcast(&testCase, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    Real dim_x = 0.0, dim_y = 0.0, dim_z = 0.0;
    Real origin_x = 0.0, origin_y = 0.0, origin_z = 0.0;
    Real extr_px = 0.0, extr_py = 0.0, extr_pz = 0.0;
    Real Re = 0.0;
    bool periodicPressureBC[3];
    DomainData *domainData = nullptr;

    switch (testCase) {
        default:
            if (IS_MAIN_PROC) cout << "Test case not listed: Abort" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            break;
        case 0:
            dim_x = 1.0;
            dim_y = 1.0;
            dim_z = 2.0;
            Re = 1000.0;

            origin_x = 0.0;
            origin_y = 0.0;
            origin_z = -1.0;

            extr_px = 0.5;
            extr_py = 0.5;
            extr_pz = 0.5;

            periodicPressureBC[0] = false;
            periodicPressureBC[1] = false;
            periodicPressureBC[2] = false;

            domainData = &testcase1DomainData;
            break;
        case 1:
            dim_x = 1.0;
            dim_y = 1.0;
            dim_z = 1.0;
            Re = 1000.0;

            origin_x = -0.5;
            origin_y = -0.5;
            origin_z = -0.5;

            extr_px = 0.0;
            extr_py = 0.0;
            extr_pz = 0.0;

            periodicPressureBC[0] = false;
            periodicPressureBC[1] = false;
            periodicPressureBC[2] = true;

            domainData = &testcase2DomainData;
        break;
    }

    enabledLogger.openSection("RUN TEST CASE");

    initSolverData(
        dim_x, dim_y, dim_z,
        origin_x, origin_y, origin_z,
        deltaT, timeSteps, Re,
        nx, ny, nz,
        npy, npz,
        periodicPressureBC,
        domainData
    );

    enabledLogger.printTitle("Data initialized");

    runSolver(extr_px, extr_py, extr_pz);

    destroySolverData();

    enabledLogger.printTitle("Data destroyed").closeSection();

    return 0;
}


int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &THIS_PROC_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &THIS_WORLD_SIZE);

    int ris = runTestCase(argc, argv);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return ris;
}
