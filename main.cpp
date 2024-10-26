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
    Real T = 0.1;
    Real deltaT = 0.0001;
    Real Re = 4700.0;

    // Define dim as side dimension of the grid (just for simplicity)
    index_t dim = 100;

    // Define number of nodes for each axis
    index_t nx = dim;
    index_t ny = dim;
    index_t nz = dim;

    // Define physical size of the problem (just for simplicity)
    Real phy_dim = 1.0;

    // Define physical size of the problem for each axis
    Real sx = phy_dim / nx;
    Real sy = phy_dim / ny;
    Real sz = phy_dim / nz;

    // Define initial velocity function
    auto initialVel = [](Real x, Real y, Real z) -> Vector {
        return {0,0,0};
    };

    // Define initial pressure function
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

        Real sdx = model.dx/2;
        Real sdy = model.dy/2;
        Real sdz = model.dz/2;

        // Iterate over face x=0 (with ghost points)
        for (index_t i = -1; i < 1; i++)
            for (index_t j = -1; j < grid.ny + 1; j++)
                for (index_t k = -1; k < grid.nz + 1; k++) {
                    Real x = static_cast<Real>(i) * model.dx;
                    Real y = static_cast<Real>(j) * model.dy;
                    Real z = static_cast<Real>(k) * model.dz;

                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx,y,z, currentTime);
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
                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx,y,z, currentTime);
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
                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx,y,z, currentTime);
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
                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx,y,z, currentTime);
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
                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx,y,z, currentTime);
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
                    grid(U, i, j, k) = ExactSolution<STANDARD>::u(x + sdx,y,z, currentTime);
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