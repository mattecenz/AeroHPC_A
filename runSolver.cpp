#include <mpi.h>
#include <cmath>
#include "Traits.hpp"
#include "utils/chronoUtils.hpp"
#include "utils/Logger.hpp"
#include "C2Decomp.hpp"
#include "L2NormCalculator.hpp"
#include "RungeKutta.hpp"
#include "DataExporter.hpp"
#include <fstream>
#include "data/SolverData.hpp"
#include "Boundaries.hpp"
#include "Initialization.hpp"
#include "vtk/VTKConverter.hpp"
#include "Interpolation.hpp"


#include "utils/printBuffer.hpp"

inline result_t runSolver(const Real extr_px, const Real extr_py, const Real extr_pz) {
    initInterpolationData(params);

    enabledLogger.printTitle("Subdomain settings")
                .printValue(5, "ID", THIS_PROC_RANK)
                .printValue(5, "IsOnTop", params.isOnTop)
                .printValue(5, "IsOnBottom", params.isOnBottom)
                .printValue(5, "IsOnLeft", params.isOnLeft)
                .printValue(5, "IsOnRight", params.isOnRight);


    enabledLogger.printTitle("Solver parameters")
                .printValue(5, "Iterations", params.timesteps)
                .printValue(5, "dT", params.dt)
                .printValue(5, "Re num", params.Re)
                .printValue(5, "Phy dim", std::to_string(params.dimX)
                                          + " x " + std::to_string(params.dimY)
                                          + " x " + std::to_string(params.dimZ));


    enabledLogger.printTitle("Grid parameters")
                .printValue(5, "nodes", std::to_string(params.loc_nX)
                                        + " x " + std::to_string(params.loc_nY)
                                        + " x " + std::to_string(params.loc_nZ));

    chrono_start(initT);
    initializeModel(rkData.model_data, 0);
    chrono_stop(initT);

    enabledLogger.printTitle("Grid initialized", initT);

    // last iteration l2Norm capture
    Real localL2NormU = 0.0;
    Real globalL2NormU = 0.0;
    Real localL2NormP = 0.0;
    Real globalL2NormP = 0.0;

    // Avg solver info
    Real avgRKTime = 0.0;
    Real avgL2Time = 0.0;
    Real avgPerf = 0.0;
    int avgCount = 0;

    // Printing variables
    index_t maxTablePrintLine = 200;
    index_t printIt = ceil(real(params.timesteps) / real(maxTablePrintLine));

    /// Start RK method ////////////////////////////////////////////////////////////////////////////////////
    enabledLogger.printTitle("Start computation")
                .openTable("Iter", {"ts", "gl2U", "gl2P", "rkT", "l2T", "TxN"});

    chrono_start(compT);

    apply_boundaries(rkData.model_data, 0, VELOCITY);

    constexpr bool ENABLE_VTK_DEBUG = true;
    const std::string vtkdir = "vtk" + std::to_string(THIS_PROC_RANK);

    if (ENABLE_VTK_DEBUG) {
        create_directories(vtkdir);
        interpolateData(rkData.model_data, 0);
        VTKConverter::exportGrid(interpData_ptr).writeFile(vtkdir + "/solution_0.vtk");
    }

    for (index_t step = 0; step < params.timesteps; ++step) {
        enabledBufferPrinter.setFile(to_string(step));

        // Compute RK step
        const Real time = real(step) * params.dt;
        chrono_start(rkTime);
        rungeKutta(time);
        chrono_stop(rkTime);

        // Prints Log every n iteration or if is the last one
        if (!((step + 1) % printIt) || step == params.timesteps - 1)
        {
            // Compute l2 norm error of velocity and pressure
            Real currentTime = time + params.dt;
            chrono_start(l2Time);
            computeL2Norm(rkData.model_data, currentTime, localL2NormU, localL2NormP);
            chrono_stop(l2Time);
            // Collect error from all processors
            globalL2NormU = 0.0;
            globalL2NormP = 0.0;
            MPI_Allreduce(&localL2NormU, &globalL2NormU, 1, Real_MPI, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&localL2NormP, &globalL2NormP, 1, Real_MPI, MPI_SUM, MPI_COMM_WORLD);
            globalL2NormU = std::sqrt(globalL2NormU);
            globalL2NormP = std::sqrt(globalL2NormP);

            // Compute per-node performance
            const Real perf = rkTime / params.grid_ndim;

            // Print log info
            enabledLogger.printTableValues(step + 1, {currentTime, globalL2NormU, globalL2NormP, rkTime, l2Time, perf});

            // Export complete vtk if needed
            if (ENABLE_VTK_DEBUG) {
                interpolateData(rkData.model_data, currentTime);
                VTKConverter::exportGrid(interpData_ptr).writeFile(vtkdir + "/solution_" + std::to_string(step) + ".vtk");
            }

            // Compute solver result info
            avgCount++;
            avgRKTime += rkTime;
            avgL2Time += l2Time;
            avgPerf += perf;
        }
    }
    chrono_stop(compT);

    enabledLogger.closeTable().printTitle("End of computation", compT);

    // Finalize solver result info
    avgRKTime /= avgCount;
    avgL2Time /= avgCount;
    avgPerf /= avgCount;

    //TODO EXPORTING
    /*
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
    */

    destroyInterpolationData();

    return {
        {"Nodes", params.glob_nX},
        {"Uerr", globalL2NormU},
        {"Perr", globalL2NormP},
        {"CompT", compT},
        {"AvgRKT", avgRKTime},
        {"AvgL2T", avgL2Time},
        {"AvgPerf", avgPerf}
    };
}
