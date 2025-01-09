#include <mpi.h>
#include <cmath>
#include "Traits.hpp"
#include "utils/chronoUtils.hpp"
#include "Logger.hpp"
#include "boundariesBuilders.cpp"
#include "C2Decomp.hpp"
#include "L2NormCalculator.hpp"
#include "RungeKutta.hpp"
#include "DataExporter.hpp"
#include <fstream>
#include "data/SolverData.hpp"

#include "printBuffer.hpp"

inline Real runSolver(const int rank, const int size,
               const boundaryDomainFunctions &boundaryDF,
               const Real extr_px, const Real extr_py, const Real extr_pz) {

    if (!rank)
        logger.openSection("Running the TestSolver")
                .printValue(5, "Iterations", params.timesteps)
                .printValue(5, "dT", params.dt)
                .printValue(5, "Re num", params.Re)
                .printValue(5, "Phy dim", std::to_string(params.dimX) + " x " + std::to_string(params.dimY) + " x " + std::to_string(params.dimZ));


    if (!rank)
        logger.printTitle("Grid created")
                .printValue(5, "nodes", std::to_string(params.loc_nX)
                                        + " x " + std::to_string(params.loc_nY)
                                        + " x " + std::to_string(params.loc_nZ));

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
    initValues(initialVel, initialPres);
    chrono_stop(initT);

    if (!rank)
        logger.printTitle("Grid initialized", initT);

    /// Define boundaries condition functions //////////////////////////////////////////////////////////////
    MPIBoundaries mpiBoundaries;
    buildMPIBoundaries(*c2d, modelStructure, mpiBoundaries, boundaryDF);

    if (!rank)
        logger.printTitle("Boundary condition set");

    /// Init variables for RK method ///////////////////////////////////////////////////////////////////////

    // last iteration l2Norm capture
    Real localL2Norm = 0.0;
    Real globalL2Norm = 0.0;

    // Performance variables
    Real perf;

    // Printing variables
    index_t maxTablePrintLine = 10;
    index_t printIt = ceil(real(params.timesteps) / real(maxTablePrintLine));

    /// Start RK method ////////////////////////////////////////////////////////////////////////////////////
    if (!rank)
        logger.printTitle("Start computation")
                .openTable("Iter", {"ts", "gl2", "rkT", "l2T", "TxN"});

    chrono_start(compT);
    mpiBoundaries.apply(model, 0);
    for (index_t step = 0; step < params.timesteps; ++step) {
        dir = "./iterations/" + to_string(step) + "/";
        c_dir();

        const Real time = real(step) * params.dt;

        chrono_start(rkTime);
        rungeKutta(time);
        chrono_stop(rkTime);

        if (!((step + 1)  % printIt) || step == params.timesteps - 1) // prints every n iteration or if is the last one
        {
            Real currentTime = time + params.dt;
            chrono_start(l2Time);
            localL2Norm = computeL2Norm(model, currentTime);
            chrono_stop(l2Time);
            perf = rkTime / params.grid_ndim;

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

    GridData interpolated_model(modelStructure);
    interpData(model, interpolated_model);

    std::array<Real, 3> originPoint = {origin_x, origin_y, origin_z};

    std::string exportFaceFilename = "solution.vtk";
    std::string exportFaceDescription = "test";
    std::vector<Real> face_points, face_vel, face_pres;
    extractFaceData(interpolated_model, face_points, face_vel, face_pres, originPoint, {0,0,0});
    writeVtkFile(exportFaceFilename, exportFaceDescription, face_points, face_vel, face_pres);

    if (!rank)
        logger.printTitle("solution written");

    std::array<Real, 3> point = {extr_px, extr_py, extr_pz};

    std::string exportLine1Filename = "profile1.dat";
    std::vector<Real> line1_points, line1_vel, line1_pres;
    extractLineData(interpolated_model, line1_points, line1_vel, line1_pres, 0, originPoint, point);
    writeDatFile(exportLine1Filename, line1_points, line1_vel, line1_pres);

    if (!rank)
        logger.printTitle("profile1 written");

    std::string exportLine2Filename = "profile2.dat";
    std::vector<Real> line2_points, line2_vel, line2_pres;
    extractLineData(interpolated_model, line2_points, line2_vel, line2_pres, 1, originPoint, point);
    writeDatFile(exportLine2Filename, line2_points, line2_vel, line2_pres);

    if (!rank)
        logger.printTitle("profile2 written");

    std::string exportLine3Filename = "profile3.dat";
    std::vector<Real> line3_points, line3_vel, line3_pres;
    extractLineData(interpolated_model, line3_points, line3_vel, line3_pres, 2, originPoint, point);
    writeDatFile(exportLine3Filename, line3_points, line3_vel, line3_pres);

    if (!rank)
        logger.printTitle("profile3 written");

    if (!rank)
        logger.closeSection().empty();

    return globalL2Norm;
}
