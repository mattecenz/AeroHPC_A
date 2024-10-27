#include <iostream>
#include <vector>
#include "L2NormCalculator.hpp"
#include "StaggeredGrid.hpp"
#include "GhostedSG.hpp"
#include "Model.hpp"
#include "BoundaryCondition.hpp"
#include "RungeKutta.hpp"
#include "VTKConverter.hpp"
#include "Chronometer.hpp"

Real testSolver(Real deltaT, index_t dim) {

    std::cout << "Running the solver" << std::endl;

    // define T & deltaT  & Re
    const Real T = 0.1;
    // const Real deltaT = 0.001;
    const Real Re = 10e6;

    // Define dim as side dimension of the grid (just for simplicity)
    // const index_t dim = 50;

    // Define number of nodes for each axis
    const index_t nx = dim;
    const index_t ny = dim;
    const index_t nz = dim;

    // Define physical size of the problem (just for simplicity)
    const Real phy_dim = 1.0;

    // Define physical size of the problem for each axis
    const Real sx = phy_dim / nx;
    const Real sy = phy_dim / ny;
    const Real sz = phy_dim / nz;

    // Define initial velocity function
    auto initialVel = [](Real x, Real y, Real z) -> Vector {
        return {0, 0, 0};
    };

    // Define initial pressure function
    // For the moment it does not work so do not care about it
    auto initialPres = [](Real x, Real y, Real z) -> Real {
        return x + y + z;
    };

    // Define spacing
    Vector spacing = {sx, sy, sz};

    // define the mesh:
    GhostedSG<STANDARD> sg({nx, ny, nz}, 1);
    cout << "Grid created" << std::endl;

    // build the model
    // initialize the grid with initial values
    Model<STANDARD> model(spacing, sg, Re, initialVel, initialPres);
    cout << "Model created, grid initialized" << endl;

    // Define test boundary condition
    BoundaryCondition<STANDARD>::Mapper testBCMapper = [&model](StaggeredGrid<STANDARD> &grid,
                                                                const TFunction &fun, Real currentTime) {

        const Real sdx = model.sdx;
        const Real sdy = model.sdy;
        const Real sdz = model.sdz;

        // Yes but when I have x=0 I need to do two things:
        // 1) set only the ghosted x points by interpolating
        // 2) set the face yz to the corresponding values

        // To understand what happens we can visualize it
        //       ---x---p---x---p---x              -> x
        // Index:  -1   0   0   1   1
        //              ^
        //           boundary

        // We want x[0] to be calculated (as that is inside our mesh so needs to be passed)
        // Because we do not know it unless we can approximate it as an int

        // We can do everything in a for loop for the x variable
        /*
        for (index_t j = -1; j < grid.ny + 1; j++){
            for (index_t k = -1; k < grid.nz + 1; k++) {
                Real x = -1                   * model.dx;
                Real y = real(j) * model.dy;
                Real z = real(k) * model.dz;

                grid(U, -1, j, k) = 2 * ExactSolution<STANDARD>::u(0., y, z, currentTime) - grid(U, 0, j,k);
                // These 2 are fine because we do not care about it
                grid(V, 0, j, k) = ExactSolution<STANDARD>::v(x, y + sdy, z, currentTime);
                grid(W, 0, j, k) = ExactSolution<STANDARD>::w(x, y, z + sdz, currentTime);
            }
        }
        // But the only tricky part is that we will need to calculate also x=0 in the rk
        */

        using Sol = ExactSolution<STANDARD>;

        // apply on face with x constant
        for (index_t j = 0; j < grid.ny; j++)
            for (index_t k = 0; k < grid.nz; k++) {
                Real y = real(j) * model.dy;
                Real z = real(k) * model.dz;

                Real x = 0;
                // On x = 0 for ghost point we have exact for U, other approximate
                grid(U, -1, j, k) = Sol::u(x, y + sdy, z + sdz, currentTime);
                grid(V, -1, j, k) = 2 * Sol::v(x, y + model.dy, z + sdz, currentTime)
                                    - grid(V, 0, j, k);
                grid(W, -1, j, k) = 2 * Sol::w(x, y + sdy, z + model.dz, currentTime)
                                    - grid(W, 0, j, k);


                x = real(grid.nx) * model.dx;
                // On x = phy_dim for domain point we have exact for U
                grid(U, grid.nx - 1, j, k) = Sol::u(x, y + sdy, z + sdz, currentTime);
                // For ghost point we have useless U, other approximate
                grid(U, grid.nx, j, k) = 0;
                grid(V, grid.nx, j, k) = 2 * Sol::v(x, y + model.dy, z + sdz, currentTime)
                                         - grid(V, grid.nx - 1, j, k);
                grid(W, grid.nx, j, k) = 2 * Sol::w(x, y + sdy, z + model.dz, currentTime)
                                         - grid(W, grid.nx - 1, j, k);
            }

        // apply on face with y constant
        for (index_t i = 0; i < grid.nx ; i++)
            for (index_t k = 0; k < grid.nz; k++) {
                Real x = real(i) * model.dx;
                Real z = real(k) * model.dz;

                Real y = 0;
                // On y = 0 for ghost point we hae exact for V, other approximate
                grid(U, i, -1, k) = 2 * Sol::u(x + model.dx, y, z + sdz, currentTime)
                                    - grid(U, i, 0, k);
                grid(V, i, -1, k) = Sol::v(x + sdx, y, z + sdz, currentTime);
                grid(W, i, -1, k) = 2 * Sol::w(x + sdx, y, z + model.dz, currentTime)
                                    - grid(W, i, 0, k);

                y = real(grid.ny) * model.dx;
                // On y = phy_dim for domain point we have exact for V
                grid(V, i, grid.ny - 1, k) = Sol::v(x + sdx, y, z + sdz, currentTime);
                // For ghost points we have useless V, other approximate
                grid(U, i, grid.ny, k) = 2 * Sol::u(x + model.dx, y, z + sdz, currentTime)
                                         - grid(U, i, grid.ny - 1, k);
                grid(V, i, grid.ny, k) = 0;
                grid(W, i, grid.ny, k) = 2 * Sol::w(x + sdx, y, z + model.dz, currentTime)
                                         - grid(W, i, grid.ny - 1, k);
            }

        // apply on face with z constant
        for (index_t i = 0; i < grid.nx; i++)
            for (index_t j = 0; j < grid.ny; j++) {
                Real x = real(i) * model.dx;
                Real y = real(j) * model.dy;

                Real z = 0;
                // On z = 0 for ghost point we have exact for W, other approximate
                grid(U, i, j, -1) = 2 * Sol::u(x + model.dx, y + sdy, z, currentTime)
                                    - grid(U, i, j, 0);
                grid(V, i, j, -1) = 2 * Sol::v(x + sdx, y + model.dy, z, currentTime)
                                    - grid(V, i, j, 0);
                grid(W, i, j, -1) = Sol::w(x + sdx, y + sdy, z, currentTime);

                z = real(grid.nz) * model.dz;
                // On z = phy_dim for domain point we have exact for W
                grid(W, i, j, grid.nz - 1) = Sol::w(x + sdx, y + sdy, z, currentTime);
                // For ghost points we gave useless W, other interpolate
                grid(U, i, j, grid.nz) = 2 * Sol::u(x + model.dx, y + sdy, z, currentTime)
                                         - grid(U, i, j, grid.nz - 1);
                grid(V, i, j, grid.nz) = 2 * Sol::v(x + sdx, y + model.dy, z, currentTime)
                                         - grid(V, i, j, grid.nz - 1);
                grid(W, i, j, grid.nz) = 0;

            }
    };

    // Define test BC function
    TFunction zero = [](Real x, Real y, Real z, Real t) {
        return 0;
    };

    // Define test BC and add it to model
    BoundaryCondition<STANDARD> inletBoundary(testBCMapper, zero);
    model.addBC(inletBoundary);
    cout << "Boundary condition set" << endl;

    // Define Buffers for RK method
    Model<STANDARD> Y2(model);
    Model<STANDARD> Y3(model);
    cout << "Buffers created" << endl;



    // Time iterations
    Real currentTime = 0.0;
    int stepCounter = 0;
    Real l2Norm = 0.0;

    model.applyBCs(currentTime);
    while (currentTime < T) {
        // call RK (obtain model at currentTime + dt)
        measure(rkTime,
                code_span(
                        rungeKutta(model, Y2, Y3, deltaT, currentTime);
                        currentTime += deltaT;
                        stepCounter++;
                )
        );

        measure(l2Time,
                code_span(
                        l2Norm = computeL2Norm<STANDARD>(model, currentTime);
                )
        );

        printf("%5d) ts %0.4f | l2 %2.7f | rkT %2.5f | l2T %2.5f\n", stepCounter, currentTime, l2Norm, rkTime, l2Time);
    }

    // Output of last iteration
    // VTKFile file = VTKConverter::exportModel(model, "testsolver output of last time iteration");
    // file.writeFile("testsolver.vtk");

    return l2Norm;
}


int main() {

    // dividing the timestep size to half
    std::vector<Real> deltaTs = {0.0001, 0.0005, 0.00025};
    std::vector<index_t> dims = {4, 8, 16, 32, 64};

    std::vector<Real> error;

    // wrt deltaT
    // for (size_t i=0; i<deltaTs.size(); ++i){
    //     Real deltaT = deltaTs[i];
    //     index_t dim = dims[0];
    //     error.push_back(testSolver(deltaT, dim));
    // }


    // wrt dim
    for (size_t i = 0; i < dims.size(); ++i) {
        Real deltaT = deltaTs[0]; // first
        index_t dim = dims[i];
        error.push_back(testSolver(deltaT, dim));
    }


    // // wrt both
    // for (size_t i=0; i<dims.size(); ++i){
    //     Real deltaT = deltaT[i];
    //     index_t dim = dims[i];
    //     error.push_back(testSolver(deltaT, dim));
    // }

    std::ofstream csvFile("output.csv");
    csvFile << "Err" << std::endl;
    for (float i: error) csvFile << i << std::endl;


    return 0;
}