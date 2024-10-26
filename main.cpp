#include <iostream>
#include "L2NormCalculator.hpp"
#include "StaggeredGrid.hpp"
#include "GhostedSG.hpp"
#include "Model.hpp"
#include "BoundaryCondition.hpp"
#include "RungeKutta.hpp"
#include "VTKConverter.hpp"

void testSolver() {

    std::cout << "Running the solver" << std::endl;

    // define T & deltaT  & Re
    const Real T = 1.;
    const Real deltaT = 0.001;
    const Real Re = 4700.0;

    // Define dim as side dimension of the grid (just for simplicity)
    const index_t dim = 50;

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
        return {0,0,0};
    };

    // Define initial pressure function
    // For the moment it does not work so do not care about it
    auto initialPres = [](Real x, Real y, Real z) -> Real {
        return x + y + z;
    };

    // Define spacing
    Vector spacing = {sx, sy, sz};

    // define the mesh:
    // instantiate from ghosted stagg grid
    // hint: second method, num of nodes in each direction, number of ghosts=1
    GhostedSG<STANDARD> sg({nx, ny, nz}, 1);
    cout << "Grid created" << std::endl;

    // build the model
    // initialize the grid with initial values
    Model<STANDARD> model(spacing, sg, Re, initialVel, initialPres);
    cout << "Model created, grid initialized" << endl;

    // Time iterations
    Real currentTime = 0.0;
    int stepCounter = 0;

    // Define test boundary condition
    BoundaryCondition<STANDARD>::Mapper testBCMapper = [&currentTime, &model](StaggeredGrid<STANDARD> &grid, const Function &fun) {

        const Real sdx = model.dx/2;
        const Real sdy = model.dy/2;
        const Real sdz = model.dz/2;

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
                Real y = static_cast<Real>(j) * model.dy;
                Real z = static_cast<Real>(k) * model.dz;

                grid(U, -1, j, k) = 2 * ExactSolution<STANDARD>::u(0., y, z, currentTime) - grid(U, 0, j,k);
                // These 2 are fine because we do not care about it
                grid(V, 0, j, k) = ExactSolution<STANDARD>::v(x, y + sdy, z, currentTime);
                grid(W, 0, j, k) = ExactSolution<STANDARD>::w(x, y, z + sdz, currentTime);
            }
        }
        // But the only tricky part is that we will need to calculate also x=0 in the rk
        */

        for (index_t i = -1; i < 1; i++)
            for (index_t j = -1; j < grid.ny + 1; j++)
                for (index_t k = -1; k < grid.nz + 1; k++) {
                    Real x = static_cast<Real>(i) * model.dx;
                    Real y = static_cast<Real>(j) * model.dy;
                    Real z = static_cast<Real>(k) * model.dz;

                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx, y, z, currentTime);
                    grid(V, i, j, k) = ExactSolution<STANDARD>::v(x, y + sdy, z, currentTime);
                    grid(W, i, j, k) = ExactSolution<STANDARD>::w(x, y, z + sdz, currentTime);
                }
        // Iterate over face x = grid.nx-1 (with ghost points)
        for (index_t i = grid.nx - 1; i < grid.nx + 1; i++)
            for (index_t j = -1; j < grid.ny + 1; j++)
                for (index_t k = -1; k < grid.nz + 1; k++) {
                    Real x = static_cast<Real>(i) * model.dx;
                    Real y = static_cast<Real>(j) * model.dy;
                    Real z = static_cast<Real>(k) * model.dz;
                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx, y, z, currentTime);
                    grid(V, i, j, k) = ExactSolution<STANDARD>::v(x, y + sdy, z, currentTime);
                    grid(W, i, j, k) = ExactSolution<STANDARD>::w(x, y, z + sdz, currentTime);
                }

        // Iterate over face y=0 (with ghost points)
        for (index_t i = -1; i < grid.nx + 1; i++)
            for (index_t j = -1; j < 1; j++)
                for (index_t k = -1; k < grid.nz + 1; k++) {
                    Real x = static_cast<Real>(i) * model.dx;
                    Real y = static_cast<Real>(j) * model.dy;
                    Real z = static_cast<Real>(k) * model.dz;
                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx, y, z, currentTime);
                    grid(V, i, j, k) = ExactSolution<STANDARD>::v(x, y + sdy, z, currentTime);
                    grid(W, i, j, k) = ExactSolution<STANDARD>::w(x, y, z + sdz, currentTime);
                }
        // Iterate over face y = grid.ny-1 (with ghost points)
        for (index_t i = -1; i < grid.nx + 1; i++)
            for (index_t j = grid.ny - 1; j < grid.ny + 1; j++)
                for (index_t k = -1; k < grid.nz + 1; k++) {
                    Real x = static_cast<Real>(i) * model.dx;
                    Real y = static_cast<Real>(j) * model.dy;
                    Real z = static_cast<Real>(k) * model.dz;
                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx, y, z, currentTime);
                    grid(V, i, j, k) = ExactSolution<STANDARD>::v(x, y + sdy, z, currentTime);
                    grid(W, i, j, k) = ExactSolution<STANDARD>::w(x, y, z + sdz, currentTime);
                }

        // Iterate over face z=0 (with ghost points)
        for (index_t i = -1; i < grid.nx + 1; i++)
            for (index_t j = -1; j < grid.ny + 1; j++)
                for (index_t k = -1; k < 1; k++) {
                    Real x = static_cast<Real>(i) * model.dx;
                    Real y = static_cast<Real>(j) * model.dy;
                    Real z = static_cast<Real>(k) * model.dz;
                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx, y, z, currentTime);
                    grid(V, i, j, k) = ExactSolution<STANDARD>::v(x, y + sdy, z, currentTime);
                    grid(W, i, j, k) = ExactSolution<STANDARD>::w(x, y, z + sdz, currentTime);
                }
        // Iterate over face z = grid.nz-1 (with ghost points)
        for (index_t i = -1; i < grid.nx + 1; i++)
            for (index_t j = -1; j < grid.ny + 1; j++)
                for (index_t k = grid.nz - 1; k < grid.nz + 1; k++) {
                    Real x = static_cast<Real>(i) * model.dx;
                    Real y = static_cast<Real>(j) * model.dy;
                    Real z = static_cast<Real>(k) * model.dz;
                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx, y, z, currentTime);
                    grid(V, i, j, k) = ExactSolution<STANDARD>::v(x, y + sdy, z, currentTime);
                    grid(W, i, j, k) = ExactSolution<STANDARD>::w(x, y, z + sdz, currentTime);
                }
    };

    // Define test BC function
    Function zero = [](Real x, Real y, Real z) {
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

    model.applyBCs(); // for T=0

    while (currentTime < T) {
        // call RK (obtain model at currentTime + dt)
        rungeKutta(model, Y2, Y3, deltaT, currentTime);

        currentTime += deltaT;
        stepCounter++;

        model.applyBCs(); // for T= currTime + dt

        /* TODO ??At this point we should have the solution on all the domain at time = currTime + dt??
             (so we can compare with exact solution at time = currTime + dt) */
        Real l2Norm = computeL2Norm<STANDARD>(model, currentTime);
        printf("%5d) t %0.4f l2 %f\n", stepCounter, currentTime, l2Norm);
    }

    // Output of last iteration
    VTKFile file = VTKConverter::exportModel(model, "testsolver output of last time iteration");
    file.writeFile("testsolver.vtk");
}

int main() {
    // try the solver
    testSolver();
    return 0;
}