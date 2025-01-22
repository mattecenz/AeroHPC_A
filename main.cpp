#include <mpi.h>
#include "Traits.hpp"
#include "DataExporter.hpp"
#include "runSolver.cpp"
#include <fstream>

using namespace std;

int runTestCase(int rank, int size, int argc, char **argv) {
    int npy, npz;
    int nx, ny, nz;
    int timeSteps;
    Real deltaT;
    int testCase;

    // Processor 0 read the file
    if (!rank) {
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

    boundaryDomainFunctions boundaryDF;
    Real extr_px, extr_py, extr_pz;

    switch (testCase) {
        default:
            if (!rank) cout << "Test case not listed: Abort" << endl;
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

            TFunction zeroFunction = [](Real x, Real y, Real z, Real t) { return real(0); };
            TFunction oneFunction = [](Real x, Real y, Real z, Real t) { return real(1); };

            boundaryDF[NORTH_FACE_ID] = {zeroFunction, zeroFunction, zeroFunction};
            boundaryDF[SOUTH_FACE_ID] = {zeroFunction, zeroFunction, zeroFunction};
            boundaryDF[EAST_FACE_ID] = {zeroFunction, zeroFunction, zeroFunction};
            boundaryDF[WEST_FACE_ID] = {zeroFunction, zeroFunction, zeroFunction};
            boundaryDF[FRONT_FACE_ID] = {zeroFunction, zeroFunction, zeroFunction};
            boundaryDF[BACK_FACE_ID] = {zeroFunction, oneFunction, zeroFunction};
            break;
    }

    runSolver(rank, size,
              npy, npz,
              nx, ny, nz,
              dim_x, dim_y, dim_z,
              deltaT, timeSteps, Re,
              boundaryDF,
              origin_x, origin_y, origin_z,
              extr_px, extr_py, extr_pz);

    return 0;
}

#include "DomainInfo.h"


int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int ris = runTestCase(rank, size, argc, argv);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return ris;
}
