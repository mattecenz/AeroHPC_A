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
    double T = 0.02;
    double deltaT = 0.0001;
    Real Re = 4700.0;
    index_t dim = 50;
    // define the mesh:
    // instantiate from ghosted stagg grid
    // hint: second method, num of nodes in each direction, number of ghosts=1
    GhostedSG<STANDARD> sg({dim, dim, dim}, 1);

    cout << "Grid created" << std::endl;

    auto initialVel = [](Real x, Real y, Real z) -> Vector {
        return {1,0,0};
    };

    auto initialPres = [](Real x, Real y, Real z) -> Real {
        return x + y + z;
    };

    Vector spacing = {0.001, 0.001, 0.001};

    // build the model
    // initialize the grid with initial values
    Model<STANDARD> model(spacing, sg, 1, initialVel, initialPres);

    cout << "Model created, grid initialized" << endl;

    BoundaryCondition<STANDARD>::Mapper inletMapping = [](StaggeredGrid<STANDARD> &grid, const Function &fun) {
        for (index_t i = -1; i < 1; i++)
            for (index_t j = -1; j < grid.ny + 1; j++)
                for (index_t k = -1; k < grid.nz + 1; k++) {
                    grid(U, k, i, j) = fun(k, i, j);
                    grid(V, k, i, j) = fun(k, i, j);
                    grid(W, k, i, j) = fun(k, i, j);
                }

        for (index_t i = grid.nx - 1; i < grid.nx + 1; i++)
            for (index_t j = -1; j < grid.ny + 1; j++)
                for (index_t k = -1; k < grid.nz + 1; k++) {
                    grid(U, k, i, j) = fun(k, i, j);
                    grid(V, k, i, j) = fun(k, i, j);
                    grid(W, k, i, j) = fun(k, i, j);
                }

        for (index_t i = -1; i < grid.nx + 1; i++)
            for (index_t j = -1; j < 1; j++)
                for (index_t k = -1; k < grid.nz + 1; k++) {
                    grid(U, k, i, j) = fun(k, i, j);
                    grid(V, k, i, j) = fun(k, i, j);
                    grid(W, k, i, j) = fun(k, i, j);
                }
        for (index_t i = -1; i < grid.nx + 1; i++)
            for (index_t j = grid.ny - 1; j < grid.ny + 1; j++)
                for (index_t k = -1; k < grid.nz + 1; k++) {
                    grid(U, k, i, j) = fun(k, i, j);
                    grid(V, k, i, j) = fun(k, i, j);
                    grid(W, k, i, j) = fun(k, i, j);
                }

        for (index_t i = -1; i < grid.nx + 1; i++)
            for (index_t j = -1; j < grid.ny + 1; j++)
                for (index_t k = -1; k < 1; k++) {
                    grid(U, k, i, j) = fun(k, i, j);
                    grid(V, k, i, j) = fun(k, i, j);
                    grid(W, k, i, j) = fun(k, i, j);
                }
        for (index_t i = -1; i < grid.nx + 1; i++)
            for (index_t j = -1; j < grid.ny + 1; j++)
                for (index_t k = grid.nz - 1; k < grid.nz + 1; k++) {
                    grid(U, k, i, j) = fun(k, i, j);
                    grid(V, k, i, j) = fun(k, i, j);
                    grid(W, k, i, j) = fun(k, i, j);
                }

    };

    Function inletFunction = [](Real x, Real y, Real z) {
        return 0;
    };

    BoundaryCondition<STANDARD> inletBoundary(inletMapping, inletFunction);
    model.addBC(inletBoundary);

    cout << "Boundary condition set" << endl;

    double currentTime = 0.0;
    int stepCounter = 0;
    Model<STANDARD> Y2(model);
    Model<STANDARD> Y3(model);


    while (currentTime < T) {
        // apply boundary cond
        model.applyBCs();
        // call RK
        rungeKutta(Re, model, Y2, Y3, deltaT, currentTime);

        Real l2Norm = computeL2Norm<STANDARD>(model, currentTime);
        printf("%d) t %0.4f: l2 %f\n", stepCounter, currentTime, l2Norm);

        currentTime += deltaT;
        stepCounter++;
    }


    VTKFile file = VTKConverter::exportModel(model, "out");
    file.writeFile("output.vtk");
}

int main() {
    // try the solver
    testSolver();
    return 0;
}