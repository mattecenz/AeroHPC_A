#include <iostream>
#include <vector>
#include "L2NormCalculator.hpp"
#include "Grid.hpp"
#include "Boundaries.hpp"
#include "RungeKutta.hpp"
#include "VTKConverter.hpp"
#include "Chronometer.hpp"

Real testSolver(Real deltaT, index_t dim) {

    // define T & deltaT  & Re
    const Real T = 1;
    // const Real deltaT = 0.001;
    const Real Re = 4700;
    // Define physical size of the problem (just for simplicity)
    const Real phy_dim = 1.0;

    std::cout << "================== Running the TestSolver ==================" << std::endl;
    cout << "  Final T: " << T << std::endl;
    cout << "       dT: " << deltaT << std::endl;
    cout << "   Re num: " << Re << endl;
    cout << "  Phy dim: " << phy_dim << " x " << phy_dim << " x " << phy_dim << std::endl;

    // Define number of nodes for each axis
    const index_t nx = dim;
    const index_t ny = dim;
    const index_t nz = dim;

    // Define physical size of the problem for each axis
    const Real sx = phy_dim / real(nx);
    const Real sy = phy_dim / real(ny);
    const Real sz = phy_dim / real(nz);

    // Define spacing
    Vector spacing = {sx, sy, sz};

    // Define nodes
    std::array<index_t,3> nodes = {nx,ny,nz};

    // define the mesh:
    Grid<STANDARD> model(nodes, spacing, 1);

    cout << "----------------------- Grid created -----------------------" << std::endl;
    cout << "  nodes: " << model.nx << " x " << model.ny << " x " << model.nz << std::endl;
    cout << " ghosts: " << model.gp << endl;


    // initialize the mesh
    // Define initial velocity function
    auto initialVel = [](Real x, Real y, Real z) -> Vector {
        return {0, 0, 0};
    };

    // Define initial pressure function
    // For the moment it does not work so do not care about it
    auto initialPres = [](Real x, Real y, Real z) -> Real {
        return x + y + z;
    };

    model.initGrid(initialVel, initialPres);
    cout << "--------------------- Grid initialized ---------------------" << std::endl;


    // Define test boundary condition
    Condition<STANDARD>::Mapper testBCMapper = [](Grid<STANDARD> &grid,
                                                  const TFunction &fun, Real currentTime) {

        const Real sdx = grid.sdx;
        const Real sdy = grid.sdy;
        const Real sdz = grid.sdz;

        using Sol = ExactSolution<STANDARD>;

        // apply on face with x constant
        for (index_t j = 0; j < grid.ny; j++)
            for (index_t k = 0; k < grid.nz; k++) {
                Real y = real(j) * grid.dy;
                Real z = real(k) * grid.dz;

                Real x = 0;
                // On x = 0 for ghost point we have exact for U, other approximate
                grid(U, -1, j, k) = Sol::u(x, y + sdy, z + sdz, currentTime);
                grid(V, -1, j, k) = 2 * Sol::v(x, y + grid.dy, z + sdz, currentTime)
                                    - grid(V, 0, j, k);
                grid(W, -1, j, k) = 2 * Sol::w(x, y + sdy, z + grid.dz, currentTime)
                                    - grid(W, 0, j, k);


                x = real(grid.nx) * grid.dx;
                // On x = phy_dim for domain point we have exact for U
                grid(U, grid.nx - 1, j, k) = Sol::u(x, y + sdy, z + sdz, currentTime);
                // For ghost point we have useless U, other approximate
                grid(U, grid.nx, j, k) = 0;
                grid(V, grid.nx, j, k) = 2 * Sol::v(x, y + grid.dy, z + sdz, currentTime)
                                         - grid(V, grid.nx - 1, j, k);
                grid(W, grid.nx, j, k) = 2 * Sol::w(x, y + sdy, z + grid.dz, currentTime)
                                         - grid(W, grid.nx - 1, j, k);
            }

        // apply on face with y constant
        for (index_t i = 0; i < grid.nx; i++)
            for (index_t k = 0; k < grid.nz; k++) {
                Real x = real(i) * grid.dx;
                Real z = real(k) * grid.dz;

                Real y = 0;
                // On y = 0 for ghost point we hae exact for V, other approximate
                grid(U, i, -1, k) = 2 * Sol::u(x + grid.dx, y, z + sdz, currentTime)
                                    - grid(U, i, 0, k);
                grid(V, i, -1, k) = Sol::v(x + sdx, y, z + sdz, currentTime);
                grid(W, i, -1, k) = 2 * Sol::w(x + sdx, y, z + grid.dz, currentTime)
                                    - grid(W, i, 0, k);

                y = real(grid.ny) * grid.dx;
                // On y = phy_dim for domain point we have exact for V
                grid(V, i, grid.ny - 1, k) = Sol::v(x + sdx, y, z + sdz, currentTime);
                // For ghost points we have useless V, other approximate
                grid(U, i, grid.ny, k) = 2 * Sol::u(x + grid.dx, y, z + sdz, currentTime)
                                         - grid(U, i, grid.ny - 1, k);
                grid(V, i, grid.ny, k) = 0;
                grid(W, i, grid.ny, k) = 2 * Sol::w(x + sdx, y, z + grid.dz, currentTime)
                                         - grid(W, i, grid.ny - 1, k);
            }

        // apply on face with z constant
        for (index_t i = 0; i < grid.nx; i++)
            for (index_t j = 0; j < grid.ny; j++) {
                Real x = real(i) * grid.dx;
                Real y = real(j) * grid.dy;

                Real z = 0;
                // On z = 0 for ghost point we have exact for W, other approximate
                grid(U, i, j, -1) = 2 * Sol::u(x + grid.dx, y + sdy, z, currentTime)
                                    - grid(U, i, j, 0);
                grid(V, i, j, -1) = 2 * Sol::v(x + sdx, y + grid.dy, z, currentTime)
                                    - grid(V, i, j, 0);
                grid(W, i, j, -1) = Sol::w(x + sdx, y + sdy, z, currentTime);

                z = real(grid.nz) * grid.dz;
                // On z = phy_dim for domain point we have exact for W
                grid(W, i, j, grid.nz - 1) = Sol::w(x + sdx, y + sdy, z, currentTime);
                // For ghost points we gave useless W, other interpolate
                grid(U, i, j, grid.nz) = 2 * Sol::u(x + grid.dx, y + sdy, z, currentTime)
                                         - grid(U, i, j, grid.nz - 1);
                grid(V, i, j, grid.nz) = 2 * Sol::v(x + sdx, y + grid.dy, z, currentTime)
                                         - grid(V, i, j, grid.nz - 1);
                grid(W, i, j, grid.nz) = 0;

            }
    };

    // Define test condition function
    TFunction zero = [](Real x, Real y, Real z, Real t) {
        return 0;
    };

    // Define test condition
    Condition<STANDARD> inletBoundary(testBCMapper, zero);

    // Define boundary conditions
    Boundaries<STANDARD> boundaries;

    // Add condition to boundaries
    boundaries.addCond(inletBoundary);
    cout << "------------------ Boundary condition set ------------------" << endl;

    // Define Buffers for RK method
    Grid<STANDARD> Y2(model.nodes, model.spacing, model.gp);
    Grid<STANDARD> Y3(model.nodes, model.spacing, model.gp);

    cout << "--------------------- Buffers created ----------------------" << endl;

    // Time
    Real currentTime = 0.0;

    // last iteration l2Norm capture
    Real l2Norm = 0.0;

    // Printing variables
    index_t iter = 0;
    index_t printIt = 100; // prints every n iterations


    boundaries.apply(model, currentTime);
    while (currentTime < T) {
        // call RK (obtain model at currentTime + dt)
        chrono_sect(rkTime,
                    code_span(
                            rungeKutta(model, Y1, Y2, Re, deltaT, currentTime, boundaries);
                            currentTime += deltaT;
                    )
        );

        chrono_sect(l2Time,
                    code_span(
                            l2Norm = computeL2Norm<STANDARD>(model, currentTime);
                    )
        );

        if (!(iter % printIt) || currentTime >= T) // prints every n iteration or if is the last one
            printf("%5ld) ts %0.4f | l2 %2.7f | rkT %2.5f | l2T %2.5f\n",
                   iter, currentTime, l2Norm, rkTime, l2Time);

        iter++;
    }

    cout << "==================== End of computation ====================" << endl << endl;
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