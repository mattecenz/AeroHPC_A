#include <mpi.h>
#include <cmath>
#include "Traits.hpp"
#include "chronoUtils.hpp"
#include "Logger.hpp"
#include "boundariesBuilders.cpp"
#include "C2Decomp.hpp"
#include "L2NormCalculator.hpp"
#include "RungeKutta.hpp"
#include "DataExporter.hpp"
#include <fstream>

using namespace std;

Real runSolver(const int npy, const int npz,
               const int nx, const int ny, const int nz,
               const Real dim_x, const Real dim_y, const Real dim_z,
               const Real deltaT, const index_t nTimeSteps, const Real Re,
               const std::string &exportFileName, const std::string &exportFileDescription) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    C2Decomp *c2d;
    bool periodicBC[3] = {true, true, true};
    c2d = new C2Decomp(nx, ny, nz, npy, npz, periodicBC);

    if (!rank)
        logger.openSection("Running the TestSolver")
                .printValue(5, "Iterations", nTimeSteps)
                .printValue(5, "dT", deltaT)
                .printValue(5, "Re num", Re)
                .printValue(5, "Phy dim", std::to_string(dim_x) + " x " + std::to_string(dim_y) + " x " + std::to_string(dim_z));

    // Define number of local nodes for each axis
    const index_t local_nx = c2d->xSize[0];
    const index_t local_ny = c2d->xSize[1];
    const index_t local_nz = c2d->xSize[2];
    Idx3 nodes = {local_nx, local_ny, local_nz};

    // Define global displacement of the grid
    const index_t px = c2d->xStart[0];
    const index_t py = c2d->xStart[1];
    const index_t pz = c2d->xStart[2];
    Idx3 displacement = {px, py, pz};

    // Define physical size of the problem for each axis
    const Real sx = dim_x / real(nx);
    const Real sy = dim_y / real(ny);
    const Real sz = dim_z / real(nz);
    Vector spacing = {sx, sy, sz};

    // Initialize the global Grid Structure
    GridStructure modelStructure(nodes, spacing, displacement, 1);

    // define the mesh:
    GridData model(modelStructure);

    if (!rank)
        logger.printTitle("Grid created")
                .printValue(5, "nodes", std::to_string(modelStructure.nx)
                                        + " x " + std::to_string(modelStructure.ny)
                                        + " x " + std::to_string(modelStructure.nz))
                .printValue(5, "ghosts", modelStructure.gp);

    /// Initialize the mesh ////////////////////////////////////////////////////////////////////////////////
    // Define initial velocity function
    auto initialVel = [](Real x, Real y, Real z) -> Vector {
        return {ExactSolution::u(x, y, z, 0), ExactSolution::v(x, y, z, 0), ExactSolution::w(x, y, z, 0)};
    };

    // Define initial pressure function
    auto initialPres = [](Real x, Real y, Real z) -> Real {
        return 0;
    };

    chrono_start(initT);
    model.initData(initialVel, initialPres);
    chrono_stop(initT);

    if (!rank)
        logger.printTitle("Grid initialized", initT);

    /// Define boundaries condition functions //////////////////////////////////////////////////////////////
    const std::vector<TFunction> boundaryFunctions = std::vector{
        ExactSolution::u,
        ExactSolution::v,
        ExactSolution::w
    };

    //MPI STUFFS

    MPIBoundaries mpiBoundaries;
    buildMPIBoundaries(*c2d, modelStructure, mpiBoundaries, boundaryFunctions);

    if (!rank)
        logger.printTitle("Boundary condition set");

    /// Create the Poisson Solver //////////////////////////////////////////////////////////////////////////
    //TODO modify parameters to a more general definition
    poissonSolver p_solver(nx, 1.0, c2d);

    if (!rank)
        logger.printTitle("Poisson solver created");

    /// Init variables for RK method ///////////////////////////////////////////////////////////////////////

    // Buffers for model data
    GridData modelBuff(modelStructure);
    // Buffers for other data
    GridStructure bufferStructure(nodes, spacing, displacement, 0);
    GridData rhsBuff(bufferStructure, false);
    GridData pressureBuff(bufferStructure, false);

    if (!rank)
        logger.printTitle("Buffers created");

    // last iteration l2Norm capture
    Real localL2Norm = 0.0;
    Real globalL2Norm = 0.0;

    // Performance variables
    Real nNodes = real(modelStructure.nx * modelStructure.ny * modelStructure.nz);
    Real perf;

    // Printing variables
    index_t printIt = 100; // prints every n iterations

    /// Start RK method ////////////////////////////////////////////////////////////////////////////////////
    if (!rank)
        logger.printTitle("Start computation")
                .openTable("Iter", {"ts", "gl2", "rkT", "l2T", "TxN"});

    chrono_start(compT);
    mpiBoundaries.apply(model, 0);
    for (index_t step = 0; step < nTimeSteps; ++step) {

        chrono_start(rkTime);
        rungeKutta(model, modelBuff, rhsBuff, pressureBuff, Re, deltaT, step, mpiBoundaries, p_solver);
        chrono_stop(rkTime);

        if (!((step + 1)  % printIt) || step == nTimeSteps - 1) // prints every n iteration or if is the last one
        {
            Real currentTime = real(step + 1) * deltaT;
            chrono_start(l2Time);
            localL2Norm = computeL2Norm(model, currentTime);
            chrono_stop(l2Time);
            perf = rkTime / nNodes;

            globalL2Norm = 0.0;
            MPI_Allreduce(&localL2Norm, &globalL2Norm, 1, Real_MPI, MPI_SUM, MPI_COMM_WORLD);
            globalL2Norm = std::sqrt(globalL2Norm);
            if (!rank)
                logger.printTableValues(step + 1, {currentTime, globalL2Norm, rkTime, l2Time, perf});
        }
    }
    chrono_stop(compT);

    if (!rank)
        logger.closeTable().printTitle("End of computation", compT);


    std::vector<Real> points, vel, pres;
    extractData(model, points, vel, pres);
    writeVtkFile(exportFileName, exportFileDescription, points, vel, pres);

    if (!rank)
        logger.printTitle("Output written");

    if (!rank)
        logger.closeSection().empty();

    return globalL2Norm;
}

void testSolver(const int npy, const int npz) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    // dividing the timestep size to half
    const std::vector<Real> deltaTs = {0.001, 0.0005, 0.00025};
    const std::vector<int> nodes = {
        /*4, */8, 16, /* 32, 64*/
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


    // wrt dim
    for (int n: nodes) {
        Real deltaT = deltaTs[0]; // first
        Real l2norm = runSolver(npy, npz,
                                n, n, n,
                                dim_x, dim_y, dim_z,
                                deltaT, timeSteps, Re,
                                "test_" + std::to_string(n) + ".vtk",
                                "Test result of grid " + std::to_string(n) + ".vtk");
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
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int npy = 2;
    int npz = 2;

    testSolver(npy, npz);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
