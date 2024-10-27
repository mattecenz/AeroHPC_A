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
    const Real Re = 4700.;

    // Define dim as side dimension of the grid (just for simplicity)
    // const index_t dim = 50;

    // Define number of nodes for each axis
    const index_t nx = dim;
    const index_t ny = dim;
    const index_t nz = dim;

    // Define physical size of the problem (just for simplicity)
    const Real phy_dim = 1.0;

    // Define physical size of the problem for each axis
    const Real sx = phy_dim / (nx-1);
    const Real sy = phy_dim / (ny-1);
    const Real sz = phy_dim / (nz-1);

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

        using Sol = ExactSolution<STANDARD>;

        // Since I am tired that nobody has a picture of the problem I will create it
        // And this will be used (with my convention on B.C), because everyone has something in mind
        // but it seems that everyone has a different idea about it

        // This is what we are using now (along one direction, the others are the same)
        //          ---x---p---x---p---x              -> x
        // Index:     -1   0   0   1   1
        //                 ^
        //              boundary
        // Spatial:  -px   0   px  dx  dx+px

        // And the boundary will be put there, if someone changes it I will personally force remove the commit
        // Then at the end we will have (Assuming we will arrive to 1)
        //          ---x---p---x---p---x              -> x
        // Index:     -1   0   0   1   1
        //                         ^
        //                      boundary
        // Spatial:          1-px  1   1+px

        // The other directions are similar
        // Which means: THE BOUNDARY IS ON THE PRESSURE. WHY? BECAUSE ITS BETTER FULL STOP.

        // Then these will be our B.C.s
        // We have to be careful, because we need to also set manually when y and z are = 0 or = grid.ny/nz
        //  Iterate over the boundary of the face x=0
        for(index_t k = 0; k < grid.nz; k++){
            const Real z = real(k) * model.dz;

            grid(U, 0, 0, k) = Sol::u(sdx, 0., z, currentTime);
            grid(V, 0, 0, k) = Sol::v(0., sdy, z, currentTime);
            grid(W, 0, 0, k) = Sol::w(0., 0., z+sdz, currentTime);

            grid(U, 0, grid.ny-1, k) = Sol::u(sdx, 1., z, currentTime);
            // The one on V is the ghost point of our face but it is on the border => useless
            grid(V, 0, grid.ny-1, k) = Sol::v(0., 1.+sdy, z, currentTime);
            grid(W, 0, grid.ny-1, k) = Sol::w(0., 1., z+sdz, currentTime);

            //NB: grid(W, 0, grid.ny-1, grid.nz-1) is useless but better to leave it
        }
        // Now iterate over the y keeping the z fixed
        // We can avoid some points as they were already calculated
        for(index_t j = 1; j < grid.ny - 1; j++){
            const Real y = real(j) * model.dy;

            grid(U, 0, j, 0) = Sol::u(sdx, j, 0., currentTime);
            grid(V, 0, j, 0) = Sol::v(0., j+sdy, 0., currentTime);
            grid(W, 0, j, 0) = Sol::w(0., j, sdz, currentTime);

            grid(U, 0, j, grid.nz-1) = Sol::u(sdx, j, 1., currentTime);
            // Here the one on V is not useless
            grid(V, 0, j, grid.nz-1) = Sol::v(0., j+sdy, 1., currentTime);
            // But the one on W is
            grid(W, 0, j, grid.nz-1) = Sol::w(0., j, 1.+sdz, currentTime);
        }
        // Do the same for y=0 and y=grid.ny-1
        // Iterate over x
        for(index_t i = 1; i < grid.nx-1; i++){
            const Real x = real(i) * model.dx;

            grid(U, i, 0, 0) = Sol::u(x+sdx, 0., 0., currentTime);
            grid(V, i, 0, 0) = Sol::v(x, sdy, 0., currentTime);
            grid(W, i, 0, 0) = Sol::w(x, 0., sdz, currentTime);

            // Do the same for z=1.
            grid(U, i, 0, grid.nz-1) = Sol::u(x+sdx, 0., 1., currentTime);
            grid(V, i, 0, grid.nz-1) = Sol::v(x, sdy, 1., currentTime);
            // The velocity on W is useless
            grid(W, i, 0, grid.nz-1) = Sol::w(x, 0., 1+sdz, currentTime);

            // Do the same now for the y=1.
            grid(U, i, grid.ny-1, 0) = Sol::u(x+sdx, 1., 0., currentTime);
            // Velocity on V is useless
            grid(V, i, grid.ny-1, 0) = Sol::v(x, 1.+sdy, 0., currentTime);
            grid(W, i, grid.ny-1, 0) = Sol::w(x, 1., sdz, currentTime);

            // Do the same for z=1.
            grid(U, i, grid.ny-1, grid.nz-1) = Sol::u(x+sdx, 1., 1., currentTime);
            // here we need only the U
            grid(V, i, grid.ny-1, grid.nz-1) = Sol::v(x, 1.+sdy, 1., currentTime);
            grid(W, i, grid.ny-1, grid.nz-1) = Sol::w(x, 1., 1.+sdz, currentTime);
        }
        // Do the same for the face at x = 1.
        for(index_t k = 0; k < grid.nz; k++){
            const Real z = real(k) * model.dz;

            // The U is not needed
            grid(U, grid.nx-1, 0, k) = Sol::u(1.+sdx, 0., z, currentTime);
            grid(V, grid.nx-1, 0, k) = Sol::v(1., sdy, z, currentTime);
            grid(W, grid.nx-1, 0, k) = Sol::w(1., 0., z+sdz, currentTime);

            grid(U, grid.nx-1, grid.ny-1, k) = Sol::u(1.+sdx, 1., z, currentTime);
            grid(V, grid.nx-1, grid.ny-1, k) = Sol::v(1., 1.+sdy, z, currentTime);
            // Here we need only w
            grid(W, grid.nx-1, grid.ny-1, k) = Sol::w(1., 1., z+sdz, currentTime);
        }
        // Now iterate over the y keeping the z fixed
        for(index_t j = 1; j < grid.ny - 1; j++){
            const Real y = real(j) * model.dy;

            // U is not needed
            grid(U, grid.nx-1, j, 0) = Sol::u(1.+sdx, j, 0., currentTime);
            grid(V, grid.nx-1, j, 0) = Sol::v(1., j+sdy, 0., currentTime);
            grid(W, grid.nx-1, j, 0) = Sol::w(1., j, sdz, currentTime);

            grid(U, grid.nx-1, j, grid.nz-1) = Sol::u(1.+sdx, j, 1., currentTime);
            // Only V is needed
            grid(V, grid.nx-1, j, grid.nz-1) = Sol::v(1., j+sdy, 1., currentTime);
            grid(W, grid.nx-1, j, grid.nz-1) = Sol::w(1., j, 1.+sdz, currentTime);
        }
        // Do inside the faces now
        for (index_t j = 1; j < grid.ny-1; j++)
            for (index_t k = 1; k < grid.nz-1; k++) {
                // This is correct
                const Real y = real(j) * model.dy;
                const Real z = real(k) * model.dz;

                // On x = 0 OUR GHOST POINT IS AT -sdx !!!
                // And to calculate that we need the exact solution in zero
                // And our y and z real coordinates do not reside on the staggered grid
                grid(U, -1, j, k) = 2*Sol::u(0., y, z, currentTime) - grid(U, 0, j, k);
                // For the x and y components they are just the exact solutions on the staggered grid
                // At x=0
                // Also in this grid we do not need the B.C on the ghost points.
                grid(V, 0, j, k) = Sol::v(0, y + sdy, z, currentTime);
                grid(W, 0, j, k) = Sol::w(0, y, z + sdz, currentTime);

                // In this way we can do also a funny thing, which is our ghost point for the x
                // will correspond to the last point of the mesh (do not believe me ? idc but you should)

                // Which means that we need to set our ghost point at just grid.nx-1

                // This is our ghost point, correctly identified but the formula is wrong
                grid(U, grid.nx - 1, j, k) = 2 * Sol::u(1., y, z, currentTime) - grid(U,grid.nx-2,j,k);
                // Our boundary also again resides in the grid.nx - 1 staggered point
                // And the ghost points are not needed in this point
                grid(V, grid.nx - 1, j, k) = Sol::v(1., y + sdy, z, currentTime);
                grid(W, grid.nx - 1, j, k) = Sol::w(1., y, z + sdz, currentTime);

                // Small note: think about in which range of x should we iterate, that should be 
                // [0, grid.nx-1)
            }
        // The exact same shit does work also on the grid with y and z 
        // (you do not believe me ? I am tired of drawing stuff so you can imagine it)
        // apply on face with y constant
        for (index_t i = 1; i < grid.nx-1; i++){
            for (index_t k = 1; k < grid.nz-1; k++) {
                const Real x = real(i) * model.dx;
                const Real z = real(k) * model.dz;

                grid(U, i, 0, k) = Sol::u(x + sdx, 0, z, currentTime);
                grid(V, i, -1, k) = 2*Sol::v(x, 0., z, currentTime) - grid(V, i, 0, k);
                grid(W, i, 0, k) = Sol::w(x, 0, z + sdz, currentTime);

                grid(U, i, grid.ny - 1, k) = Sol::u(x + sdx, 1., z, currentTime);
                grid(V, i, grid.ny - 1, k) = 2 * Sol::v(x, 1., z, currentTime) - grid(V,i,grid.ny-2,k);
                grid(W, i, grid.ny - 1, k) = Sol::w(x, 1., z + sdz, currentTime);
            }
        }
        // apply on face with z constant
        for (index_t i = 1; i < grid.nx-1; i++)
            for (index_t j = 1; j < grid.ny-1; j++) {
                const Real x = real(i) * model.dx;
                const Real y = real(j) * model.dy;

                grid(U, i, j, 0) = Sol::u(x + sdx, y, 0., currentTime);
                grid(V, i, j, 0) = Sol::v(x, y + sdy, 0., currentTime);
                grid(W, i, j, -1) = 2*Sol::w(x, y, 0., currentTime) - grid(W, i, j, 0);

                grid(U, i, j, grid.nz - 1) = Sol::u(x + sdx, y, 1., currentTime);
                grid(V, i, j, grid.nz - 1) = Sol::v(x, y + sdy, 1., currentTime);
                grid(W, i, j, grid.nz - 1) = 2 * Sol::w(x, y, 1., currentTime) - grid(W,i,j,grid.nz-2);
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

        // if(stepCounter==979){
        //     std::cout<<std::endl<<"Printing model at timestep: "<<stepCounter<<std::endl<<std::endl;
        //     // Print the whole matrix why not
        //     for(index_t i=-1;i<model.grid.nx+1;++i){
        //         for(index_t j=-1;j<model.grid.ny+1;++j){
        //             for(index_t k=-1;k<model.grid.nz+1;++k){
        //                 std::cout<<"("<<U<<","<<i<<","<<j<<","<<k<<"): "<<model.grid(U,i,j,k)<<std::endl;
        //                 std::cout<<"("<<V<<","<<i<<","<<j<<","<<k<<"): "<<model.grid(V,i,j,k)<<std::endl;
        //                 std::cout<<"("<<W<<","<<i<<","<<j<<","<<k<<"): "<<model.grid(W,i,j,k)<<std::endl;
        //             }
        //         }
        //     }
        // }

        printf("%5d) ts %0.4f | l2 %2.7f | rkT %2.5f | l2T %2.5f\n", stepCounter, currentTime, l2Norm, rkTime, l2Time);
    }

    // Output of last iteration
    // VTKFile file = VTKConverter::exportModel(model, "testsolver output of last time iteration");
    // file.writeFile("testsolver.vtk");

    return l2Norm;
}


int main() {

    // dividing the timestep size to half
    std::vector<Real> deltaTs = {0.0001, 0.0005, 0.00025, 0.000125};
    std::vector<index_t> dims = {4, 8, 16, 32};

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


    // wrt both
    // for (size_t i=0; i<dims.size(); ++i){
    //     Real deltaT = deltaTs[i];
    //     index_t dim = dims[i];
    //     error.push_back(testSolver(deltaT, dim));
    // }

    std::ofstream csvFile("output.csv");
    csvFile << "Err" << std::endl;
    for (float i: error) csvFile << i << std::endl;


    return 0;
}