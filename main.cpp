#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <fstream>
#include <assert.h>

#include "Traits.hpp"
#include "VTKConverter.hpp"
#include "chronoUtils.hpp"
#include "Logger.hpp"
#include "boundariesBuilders.cpp"
#include "C2Decomp.hpp"
#include "L2NormCalculator.hpp"
#include "RungeKutta.hpp"

using namespace std;


Real testSolver(Real deltaT, index_t dim) {

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    C2Decomp *c2d;
    int pRow = 1, pCol = size;
    bool periodicBC[3] = {true, true, true};
    c2d = new C2Decomp(dim, dim, dim, pRow, pCol, periodicBC);

    // define T & deltaT  & Re
    constexpr Real T = 0.1;
    constexpr Real Re = 1000;
    // Define physical size of the problem (just for simplicity)
    constexpr Real phy_dim = 1.0;

    if (!rank)
        logger.openSection("Running the TestSolver")
            .printValue(5, "Final T", T)
            .printValue(5, "dT", deltaT)
            .printValue(5, "Re num", Re)
            .printValue(5, "Phy dim", std::to_string(phy_dim) + " x " + std::to_string(phy_dim) + " x " + std::to_string(phy_dim));

    // Define number of nodes for each axis
    const index_t nx = c2d->xSize[0];
    const index_t ny = c2d->xSize[1];
    const index_t nz = c2d->xSize[2];
    Idx3 nodes = {nx, ny, nz};

    // Define global displacement of the grid
    // NOT SURE
    const index_t px = c2d->xStart[0];
    const index_t py = c2d->xStart[1];
    const index_t pz = c2d->xStart[2];
    Idx3 displacement = {px, py, pz};

    // Define physical size of the problem for each axis
    const Real sx = phy_dim / real(nx);
    const Real sy = phy_dim / real(ny);
    const Real sz = phy_dim / real(nz);
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
    // For the moment it does not work so do not care about it
    auto initialPres = [](Real x, Real y, Real z) -> Real {
        return x + y + z;
    };

    chrono_start(initT);
    model.initData(initialVel, initialPres);
    chrono_stop(initT);

    if (!rank)
        logger.printTitle("Grid initialized", initT);

    /// Define boundaries condition functions //////////////////////////////////////////////////////////////
    const std::vector<TFunction> boundaryFunctions = std::vector{ExactSolution::u,
                                                                 ExactSolution::v,
                                                                 ExactSolution::w};

    //MPI STUFFS
    
    MPIBoundaries mpiBoundaries;
    buildMPIBoundaries(*c2d, modelStructure, mpiBoundaries, boundaryFunctions);
    

    if (!rank)
        logger.printTitle("Boundary condition set");

    /// Init variables for RK method ///////////////////////////////////////////////////////////////////////

    // Buffers
    GridData modelBuff(modelStructure);
    GridData rhsBuff(modelStructure);
    if (!rank)
        logger.printTitle("Buffers created");

    // Time
    Real currentTime = 0.0;

    // last iteration l2Norm capture
    Real localL2Norm = 0.0;
    Real globalL2Norm = 0.0;

    // Performance variables
    Real nNodes = real(modelStructure.nx * modelStructure.ny * modelStructure.nz);
    Real perf;

    // Printing variables
    index_t iter = 0;
    index_t printIt = 100; // prints every n iterations

    /// Start RK method ////////////////////////////////////////////////////////////////////////////////////
    if (!rank)
        logger.printTitle("Start computation")
            .openTable("Iter", {"ts", "l2", "rkT", "l2T", "TxN"});

    chrono_start(compT);
    mpiBoundaries.apply(model, currentTime);
    while (currentTime < T) {
        // call RK (obtain model at currentTime + dt)
        chrono_start(rkTime);
        rungeKutta(model, modelBuff, rhsBuff, Re, deltaT, currentTime, mpiBoundaries);
        currentTime += deltaT;
        chrono_stop(rkTime);

        if (!(iter % printIt) || currentTime >= T) // prints every n iteration or if is the last one
        {
            chrono_start(l2Time);
            localL2Norm = computeL2Norm(model, currentTime);
            chrono_stop(l2Time);
            perf = rkTime / nNodes;

            MPI_Allreduce(&localL2Norm, &globalL2Norm, 1, Real_MPI, MPI_SUM, MPI_COMM_WORLD);
            globalL2Norm = std::sqrt(globalL2Norm);
            if (!rank)
                logger.printTableValues(iter, {currentTime, localL2Norm, rkTime, l2Time, perf});
        }
        iter++;
    }
    chrono_stop(compT);

    if (!rank)
        logger.closeTable()
            .printTitle("End of computation", compT)
            .closeSection()
            .empty();

    return globalL2Norm;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    // dividing the timestep size to half
    std::vector<Real> deltaTs = {0.0001, 0.0005, 0.00025};
    std::vector<index_t> dims = {4, 8, 16, 32, 64};
    //std::vector<index_t> dims = {8};

    std::vector<Real> error;

    // wrt deltaT
    // for (size_t i=0; i<deltaTs.size(); ++i){
    //     Real deltaT = deltaTs[i];
    //     index_t dim = dims[3];
    //     error.push_back(testSolver(deltaT, dim));
    // }


    // wrt dim
    for (long dim: dims) {
        Real deltaT = deltaTs[0]; // first
        error.push_back(testSolver(deltaT, dim));
    }


    // // wrt both
    // for (size_t i=0; i<dims.size(); ++i){
    //     Real deltaT = deltaT[i];
    //     index_t dim = dims[i];
    //     error.push_back(testSolver(deltaT, dim));
    // }

    std::ofstream csvFile("output.csv");
    csvFile << "step,error" << std::endl;
    for (int i = 0; i < dims.size(); ++i) csvFile << dims[i] << "," << error[i] << std::endl;

    MPI_Finalize();

    return 0;
}