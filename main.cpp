#include <vector>
#include "L2NormCalculator.hpp"
#include "GridData.hpp"
#include "Boundaries.hpp"
#include "Conditions.hpp"
#include "RungeKutta.hpp"
#include "VTKConverter.hpp"
#include "Chronometer.hpp"
#include "Logger.hpp"
#include "GridStructure.hpp"


Real testSolver(Real deltaT, index_t dim) {
    // define T & deltaT  & Re
    constexpr Real T = 1;
    constexpr Real Re = 4000;
    // Define physical size of the problem (just for simplicity)
    constexpr Real phy_dim = 1.0;

    logger.openSection("Running the TestSolver");
    logger.printValue(5, "Final T", T);
    logger.printValue(5, "dT", deltaT);
    logger.printValue(5, "Re num", Re);
    logger.printValue(5, "Phy dim", to_string(phy_dim) + " x " + to_string(phy_dim) + " x " + to_string(phy_dim));

    // Define number of nodes for each axis
    const index_t nx = dim;
    const index_t ny = dim;
    const index_t nz = dim;
    Idx3 nodes = {nx,ny,nz};

    // Define global displacement of the grid
    const index_t gx = 0;
    const index_t gy = 0;
    const index_t gz = 0;
    Idx3 displacement = {gx,gy,gz};

    // Define physical size of the problem for each axis
    const Real sx = phy_dim / real(nx);
    const Real sy = phy_dim / real(ny);
    const Real sz = phy_dim / real(nz);
    Vector spacing = {sx, sy, sz};

    // Initialize the global Grid Structure
    GridStructure modelStructure(nodes, spacing, displacement, 1);

    // define the mesh:
    GridData model(modelStructure);

    logger.printTitle("Grid created");
    logger.printValue(5, "nodes", to_string(modelStructure.nx) + " x " + to_string(modelStructure.ny) + " x " + to_string(modelStructure.nz));
    logger.printValue(5, "ghosts", modelStructure.gp);

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

    logger.printTitle("Grid initialized", initT);

    /// Define boundaries condition functions //////////////////////////////////////////////////////////////
    const std::vector<TFunction> boundaryFunctions = std::vector{ExactSolution::u,
                                                                 ExactSolution::v,
                                                                 ExactSolution::w};

    /// Define face mappers ////////////////////////////////////////////////////////////////////////////////
    // North
    PhysicalCondition::Mapper northFace = [](GridData &grid,
                                             const Real currentTime,
                                             const std::vector<TFunction> &functions) {
        const Real sdx = grid.structure.sdx;
        const Real sdy = grid.structure.sdy;
        const Real sdz = grid.structure.sdz;

        const TFunction exactU = functions[0];
        const TFunction exactV = functions[1];
        const TFunction exactW = functions[2];

        const index_t j = 0;
        const Real y = real(j) * grid.structure.dy;

        // apply on face with y constant

        for (index_t k = 0; k < grid.structure.nz; k++)
            for (index_t i = 0; i < grid.structure.nx; i++) {
                Real x = real(i) * grid.structure.dx;
                Real z = real(k) * grid.structure.dz;

                // On y = 0 for ghost point we hae exact for V, other approximate
                grid.U(i, j - 1, k) = 2 * exactU(x + grid.structure.dx, y, z + sdz, currentTime)
                                      - grid.U(i, j, k);
                grid.V(i, j - 1, k) = exactV(x + sdx, y, z + sdz, currentTime);
                grid.W(i, j - 1, k) = 2 * exactW(x + sdx, y, z + grid.structure.dz, currentTime)
                                      - grid.W(i, j, k);
            }


    };
    // South
    PhysicalCondition::Mapper southFace = [](GridData &grid,
                                             const Real currentTime,
                                             const std::vector<TFunction> &functions) {

        const Real sdx = grid.structure.sdx;
        const Real sdy = grid.structure.sdy;
        const Real sdz = grid.structure.sdz;

        const TFunction exactU = functions[0];
        const TFunction exactV = functions[1];
        const TFunction exactW = functions[2];

        const index_t j = grid.structure.ny;
        const Real y = real(j) * grid.structure.dy;

        // apply on face with y constant

        for (index_t k = 0; k < grid.structure.nz; k++)
            for (index_t i = 0; i < grid.structure.nx; i++) {
                Real x = real(i) * grid.structure.dx;
                Real z = real(k) * grid.structure.dz;

                // On y = phy_dim for domain point we have exact for V
                grid.V(i, j - 1, k) = exactV(x + sdx, y, z + sdz, currentTime);
                // For ghost points we have useless V, other approximate

                grid.U(i, j, k) = 2 * exactU(x + grid.structure.dx, y, z + sdz, currentTime)
                                  - grid.U(i, j - 1, k);
                grid.V(i, j, k) = 0;
                grid.W(i, j, k) = 2 * exactW(x + sdx, y, z + grid.structure.dz, currentTime)
                                  - grid.W(i, j - 1, k);
            }
    };

    // East
    PhysicalCondition::Mapper eastFace = [](GridData &grid,
                                            const Real currentTime,
                                            const std::vector<TFunction> &functions) {
        const Real sdx = grid.structure.sdx;
        const Real sdy = grid.structure.sdy;
        const Real sdz = grid.structure.sdz;

        const TFunction exactU = functions[0];
        const TFunction exactV = functions[1];
        const TFunction exactW = functions[2];

        const index_t k = grid.structure.nz;
        const Real z = real(k) * grid.structure.dz;

        for (index_t j = 0; j < grid.structure.ny; j++)
            for (index_t i = 0; i < grid.structure.nx; i++) {
                Real x = real(i) * grid.structure.dx;
                Real y = real(j) * grid.structure.dy;

                // On z = phy_dim for domain point we have exact for W
                grid.W(i, j, k - 1) = exactW(x + sdx, y + sdy, z, currentTime);
                // For ghost points we gave useless W, other interpolate
                grid.U(i, j, k) = 2 * exactU(x + grid.structure.dx, y + sdy, z, currentTime)
                                  - grid.U(i, j, k - 1);
                grid.V(i, j, k) = 2 * exactV(x + sdx, y + grid.structure.dy, z, currentTime)
                                  - grid.V(i, j, k - 1);
                grid.W(i, j, k) = 0;

            }

    };
    // West
    PhysicalCondition::Mapper westFace = [](GridData &grid,
                                            const Real currentTime,
                                            const std::vector<TFunction> &functions) {
        const Real sdx = grid.structure.sdx;
        const Real sdy = grid.structure.sdy;
        const Real sdz = grid.structure.sdz;

        const TFunction exactU = functions[0];
        const TFunction exactV = functions[1];
        const TFunction exactW = functions[2];

        const index_t k = 0;
        const Real z = real(k) * grid.structure.dz;

        for (index_t j = 0; j < grid.structure.ny; j++)
            for (index_t i = 0; i < grid.structure.nx; i++) {
                Real x = real(i) * grid.structure.dx;
                Real y = real(j) * grid.structure.dy;

                // On z = 0 for ghost point we have exact for W, other approximate
                grid.U(i, j, k - 1) = 2 * exactU(x + grid.structure.dx, y + sdy, z, currentTime)
                                      - grid.U(i, j, k);
                grid.V(i, j, k - 1) = 2 * exactV(x + sdx, y + grid.structure.dy, z, currentTime)
                                      - grid.V(i, j, k);
                grid.W(i, j, k - 1) = exactW(x + sdx, y + sdy, z, currentTime);
            }

    };

    // Front
    PhysicalCondition::Mapper frontFace = [](GridData &grid,
                                             const Real currentTime,
                                             const std::vector<TFunction> &functions) {
        const Real sdx = grid.structure.sdx;
        const Real sdy = grid.structure.sdy;
        const Real sdz = grid.structure.sdz;

        const TFunction exactU = functions[0];
        const TFunction exactV = functions[1];
        const TFunction exactW = functions[2];

        const index_t i = 0;
        const Real x = real(i) * grid.structure.dx;

        for (index_t k = 0; k < grid.structure.nz; k++)
            for (index_t j = 0; j < grid.structure.ny; j++) {
                Real y = real(j) * grid.structure.dy;
                Real z = real(k) * grid.structure.dz;

                // On x = 0 for ghost point we have exact for U, other approximate
                grid.U(i - 1, j, k) = exactU(x, y + sdy, z + sdz, currentTime);
                grid.V(i - 1, j, k) = 2 * exactV(x, y + grid.structure.dy, z + sdz, currentTime)
                                      - grid.V(i, j, k);
                grid.W(i - 1, j, k) = 2 * exactW(x, y + sdy, z + grid.structure.dz, currentTime)
                                      - grid.W(i, j, k);
            }

    };

    // Back
    PhysicalCondition::Mapper backFace = [](GridData &grid,
                                            const Real currentTime,
                                            const std::vector<TFunction> &functions) {
        const Real sdx = grid.structure.sdx;
        const Real sdy = grid.structure.sdy;
        const Real sdz = grid.structure.sdz;

        const TFunction exactU = functions[0];
        const TFunction exactV = functions[1];
        const TFunction exactW = functions[2];

        const index_t i = grid.structure.nx;
        const Real x = real(i) * grid.structure.dx;

        for (index_t k = 0; k < grid.structure.nz; k++)
            for (index_t j = 0; j < grid.structure.ny; j++) {
                Real y = real(j) * grid.structure.dy;
                Real z = real(k) * grid.structure.dz;

                // On x = phy_dim for domain point we have exact for U
                grid.U(i - 1, j, k) = exactU(x, y + sdy, z + sdz, currentTime);
                // For ghost point we have useless U, other approximate
                grid.U(i, j, k) = 0;
                grid.V(i, j, k) = 2 * exactV(x, y + grid.structure.dy, z + sdz, currentTime)
                                  - grid.V(i - 1, j, k);
                grid.W(i, j, k) = 2 * exactW(x, y + sdy, z + grid.structure.dz, currentTime)
                                  - grid.W(i - 1, j, k);
            }
    };

    /// Define boundary conditions /////////////////////////////////////////////////////////////////////////
    PhysicalCondition northCond(northFace, boundaryFunctions);
    PhysicalCondition southCond(southFace, boundaryFunctions);
    PhysicalCondition eastCond(eastFace, boundaryFunctions);
    PhysicalCondition westCond(westFace, boundaryFunctions);
    PhysicalCondition frontCond(frontFace, boundaryFunctions);
    PhysicalCondition backCond(backFace, boundaryFunctions);

    /// Add boundary conditions to collection //////////////////////////////////////////////////////////////
    Boundaries boundaries({
                                  &northCond,
                                  &southCond,
                                  &eastCond,
                                  &westCond,
                                  &frontCond,
                                  &backCond
                          });

    logger.printTitle("Boundary condition set");

    /// Init variables for RK method ///////////////////////////////////////////////////////////////////////

    // Buffers
    GridData modelBuff(modelStructure);
    GridData rhsBuff(modelStructure);
    logger.printTitle("Buffers created");

    // Time
    Real currentTime = 0.0;

    // last iteration l2Norm capture
    Real l2Norm = 0.0;

    // Performance variables
    Real nNodes = real(modelStructure.nx * modelStructure.ny * modelStructure.nz);
    Real perf;

    // Printing variables
    index_t iter = 0;
    index_t printIt = 100; // prints every n iterations

    /// Start RK method ////////////////////////////////////////////////////////////////////////////////////
    logger.printTitle("Start computation");
    logger.openTable("Iter", {"ts", "l2", "rkT", "l2T", "TxN"});

    chrono_start(compT);
    boundaries.apply(model, currentTime);
    while (currentTime < T) {
        // call RK (obtain model at currentTime + dt)
        chrono_start(rkTime);
        rungeKutta(model, modelBuff, rhsBuff, Re, deltaT, currentTime, boundaries);
        currentTime += deltaT;
        chrono_stop(rkTime);


        if (!(iter % printIt) || currentTime >= T) // prints every n iteration or if is the last one
        {
            chrono_start(l2Time);
            l2Norm = computeL2Norm(model, currentTime);
            chrono_stop(l2Time);
            perf = rkTime / nNodes;
            logger.printTableValues(iter, {currentTime, l2Norm, rkTime, l2Time, perf});
        }
        iter++;
    }
    chrono_stop(compT);

    logger.closeTable();
    logger.printTitle("End of computation", compT);
    logger.closeSection();
    logger.empty();

    return l2Norm;
}


int main() {
    // dividing the timestep size to half
    std::vector<Real> deltaTs = {0.001, 0.0005, 0.00025};
    std::vector<index_t> dims = {4, 8, 16, 32, 64};

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

    return 0;
}