#include <iostream>
#include <vector>

#include "L2NormCalculator.hpp"
#include "StaggeredGrid.hpp"
#include "GhostedSG.hpp"
#include "Model.hpp"
#include "BoundaryCondition.hpp"
#include "RungeKutta.hpp"
#include "VTKConverter.hpp"

Real testSolver(Real deltaT, index_t dim) {

    std::cout << "Running the solver" << std::endl;

    // define T & deltaT  & Re
    const Real T = 0.1;
    // const Real deltaT = 0.001;
    const Real Re = 4700.0;

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

    Real l2Norm = 0.0;
    while (currentTime < T) {
        // call RK (obtain model at currentTime + dt)
        rungeKutta(model, Y2, Y3, deltaT, currentTime);

        currentTime += deltaT;
        stepCounter++;

        model.applyBCs(); // for T= currTime + dt

        if(currentTime >= T){ 
            l2Norm = computeL2Norm<STANDARD>(model, currentTime);
            printf("%5d) t %0.4f l2 %f\n", stepCounter, currentTime, l2Norm);
        }
        // l2Norm = computeL2Norm<STANDARD>(model, currentTime);
        // printf("%5d) t %0.4f l2 %f\n", stepCounter, currentTime, l2Norm);
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
    for (size_t i=0; i<dims.size(); ++i){
        Real deltaT = deltaTs[0]; // first
        index_t dim = dims[i];
        testSolver(deltaT, dim);
    }


    // // wrt both
    // for (size_t i=0; i<dims.size(); ++i){
    //     Real deltaT = deltaT[i];
    //     index_t dim = dims[i];
    //     testSolver(deltaT, dim);
    // }

    for (size_t i=0; i<error.size(); ++i) std::cout << "err: " << error[i] << std::endl; 
    std::ofstream csvFile("output.csv");

    return 0;
}