#include <iostream>
#include "L2NormCalculator.hpp"
#include "StaggeredGrid.hpp"
#include "GhostedSG.hpp"
#include "Model.hpp"
#include "BoundaryCondition.hpp"
#include "RungeKutta.hpp"

// Declare the RungeKutta function if not already declared in the included headers
void RungeKutta(int order, Real Re, Model<Addressing_T::STANDARD> &model, Model<Addressing_T::STANDARD> &model_out, double deltaT, double currentTime);

void solver() {

    std::cout << "Running the solver" << std::endl;

    // define T & deltaT  & Re
    double T = 0.1;
    double deltaT = 0.0001;
    Real Re = 47000.0;

    // define the mesh:
    // instantiate from ghosted stagg grid
    // hint: second method, num of nodes in each direction, number of ghosts=1
    GhostedSG<Addressing_T::STANDARD> sg({5, 5, 5}, 1);

    auto initialVel = [](Real x, Real y, Real z) -> Vector {
        return {x, y, z};
    };

    auto initialPres = [](Real x, Real y, Real z) -> Real {
        return x + y + z;
    };

    Vector spacing = {1, 1, 1};

    // build the model
    // initialize the grid with initial values
    Model<Addressing_T::STANDARD> model(spacing, sg, 1, initialVel, initialPres);

    BoundaryCondition<Addressing_T::STANDARD>::Mapper inletMapping = [](StaggeredGrid<Addressing_T::STANDARD> &grid, const Function &fun) {
        for (int i = 0; i < grid.ny; i++)
            for (int j = 0; j < grid.nz; j++) {
                grid(U, 0, i, j) = fun(0, i, j);
                grid(V, 0, i, j) = fun(0, i, j);
                grid(W, 0, i, j) = fun(0, i, j);
            }
    };

    Function inletFunction = [](Real x, Real y, Real z) {
        return 0;
    };

    BoundaryCondition<Addressing_T::STANDARD> inletBoundary(inletMapping, inletFunction);
    model.addBC(inletBoundary);

    // ERROR
    model.applyBCs();


    // couldnt test, try later
//     double time = 0.001;
//     double l2Norm = computeL2Norm<Addressing_T::STANDARD>(model, time);

    //Output Model:
    Model<Addressing_T::STANDARD> model_out(model);

    double currentTime = 0.0;
    int stepCounter = 0;
    while (currentTime < T) {
        // apply boundary cond
        // call RK
    RungeKutta(5, Re, model, model_out, deltaT, currentTime);
        currentTime += deltaT;
        stepCounter++;
    }
}

int main() {

    // try the solver
    solver();
    return 0;
}