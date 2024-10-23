#include <iostream>
#include "StaggeredGrid.hpp"
#include "GhostedSG.hpp"
#include "Model.hpp"

void solver() {

    std::cout << "Runing the solver" << std::endl;
    
    // define T & deltaT
    double T = 0.1;
    double deltaT = 0.0001;

    // define the mesh:
    // instantiate from ghosted stagg grid
    // hint: second method, num of nodes in each direction, number of ghosts=1
    GhostedSG<Addressing_T::STANDARD> sg({5, 5, 5},1);

    auto initialVel = [](Real x, Real y, Real z) -> Vector {
        return {x,y,z};
    };

    auto initialPres = [](Real x, Real y, Real z) -> Real {
        return x+y+z;
    };

    Vector spacing = {1,1,1};

    // build the model
    // initialize the grid with initial values
    Model<Addressing_T::STANDARD> model(spacing, sg, 1, initialVel, initialPres);

    auto inletMapping = [](StaggeredGrid<Addressing_T::STANDARD> &grid, Function &fun){
        for(int i = 0; i<grid.ny; i++)
            for(int j= 0; i<grid.nz; j++) {
                grid(U, 0, i, j) = fun(0,i,j);
                grid(V, 0, i, j) = fun(0,i,j);
                grid(W, 0, i, j) = fun(0,i,j);
            }
    };

    auto inletFunction = [](Real x, Real y, Real z) -> Real {
        return 0;
    };

    BoundaryCondition<Addressing_T::STANDARD> inletBoundary(inletMapping, inletFunction);
    model.addBC(inletBoundary);

    // causing errors
    // model.applyBCs();


    double currentTime = 0.0;
    unsigned int stepCounter = 0;
    while (currentTime < T){
        // apply boundary cond
        // call RK

        currentTime += deltaT;
        stepCounter ++;
    }
}

int main() {
    
    solver();
    return 0;
}