#include <iostream>
#include <vector>
#include "L2NormCalculator.hpp"
#include "Grid.hpp"
#include "Boundaries.hpp"
#include "VTKConverter.hpp"
#include "Chronometer.hpp"
#include "Operators.hpp"





class TestCaseIdentity : public VectorialFunction
{
public:
    Real u(Real x, Real y, Real z) override
    {
        return std::sin(x) * std::cos(y) * std::sin(z);
    }

    Real v(Real x, Real y, Real z) override
    {
        return std::cos(x) * std::sin(y) * std::sin(z);
    }

    Real w(Real x, Real y, Real z) override
    {
        return 2 * std::cos(x) * std::cos(y) * std::cos(z);
    }
};


class TestCaseLaplacian : public VectorialFunction
{
public:
    Real u(Real x, Real y, Real z) override
    {
        return - 3 * std::sin(x) * std::cos(y) * std::sin(z);
    }

    Real v(Real x, Real y, Real z) override
    {
        return -3 * std::cos(x) * std::sin(y) * std::sin(z);
    }

    Real w(Real x, Real y, Real z) override
    {
        return -6 * std::cos(x) * std::cos(y) * std::cos(z);
    }
};


class TestCaseConvection : public VectorialFunction
{
public:
    Real u(Real x, Real y, Real z) override
    {
        Real f1 = std::sin(x) * std::cos(y) * std::sin(z);
        Real f2 = std::cos(x) * std::sin(y) * std::sin(z);
        Real f3 = 2 * std::cos(x) * std::cos(y) * std::cos(z);

        Real df1_dx = std::cos(x) * std::cos(y) * std::sin(z);
        Real df1_dy = -std::sin(x) * std::sin(y) * std::sin(z);
        Real df1_dz = std::sin(x) * std::cos(y) * std::cos(z);

        return f1 * df1_dx + f2 * df1_dy + f3 * df1_dz;
    }

    Real v(Real x, Real y, Real z) override
    {
        Real f1 = std::sin(x) * std::cos(y) * std::sin(z);
        Real f2 = std::cos(x) * std::sin(y) * std::sin(z);
        Real f3 = 2 * std::cos(x) * std::cos(y) * std::cos(z);

        Real df2_dx = -std::sin(x) * std::sin(y) * std::sin(z);
        Real df2_dy = std::cos(x) * std::cos(y) * std::sin(z);
        Real df2_dz = std::cos(x) * std::sin(y) * std::cos(z);

        return f1 * df2_dx + f2 * df2_dy + f3 * df2_dz;
    }

    Real w(Real x, Real y, Real z) override
    {
        Real f1 = std::sin(x) * std::cos(y) * std::sin(z);
        Real f2 = std::cos(x) * std::sin(y) * std::sin(z);
        Real f3 = 2 * std::cos(x) * std::cos(y) * std::cos(z);

        Real df3_dx = -2 * std::sin(x) * std::cos(y) * std::cos(z);
        Real df3_dy = -2 * std::cos(x) * std::sin(y) * std::cos(z);
        Real df3_dz = -2 * std::cos(x) * std::cos(y) * std::sin(z);

        return f1 * df3_dx + f2 * df3_dy + f3 * df3_dz;
    }
};


class TestCaseComposed : public VectorialFunction
{
private:
    const Real Re;

public: 
    TestCaseComposed(Real Re_) : Re(Re_) {}


    Real u(Real x, Real y, Real z) override
    {
        Real f1 = std::sin(x) * std::cos(y) * std::sin(z);
        Real f2 = std::cos(x) * std::sin(y) * std::sin(z);
        Real f3 = 2 * std::cos(x) * std::cos(y) * std::cos(z);

        Real df1_dx = std::cos(x) * std::cos(y) * std::sin(z);
        Real df1_dy = -std::sin(x) * std::sin(y) * std::sin(z);
        Real df1_dz = std::sin(x) * std::cos(y) * std::cos(z);

        return f1 * df1_dx + f2 * df1_dy + f3 * df1_dz - (3 / Re) * f1;
    }

    Real v(Real x, Real y, Real z) override
    {
        Real f1 = std::sin(x) * std::cos(y) * std::sin(z);
        Real f2 = std::cos(x) * std::sin(y) * std::sin(z);
        Real f3 = 2 * std::cos(x) * std::cos(y) * std::cos(z);

        Real df2_dx = -std::sin(x) * std::sin(y) * std::sin(z);
        Real df2_dy = std::cos(x) * std::cos(y) * std::sin(z);
        Real df2_dz = std::cos(x) * std::sin(y) * std::cos(z);

        return f1 * df2_dx + f2 * df2_dy + f3 * df2_dz - (3 / Re) * f2;
    }

    Real w(Real x, Real y, Real z) override
    {
        Real f1 = std::sin(x) * std::cos(y) * std::sin(z);
        Real f2 = std::cos(x) * std::sin(y) * std::sin(z);
        Real f3 = 2 * std::cos(x) * std::cos(y) * std::cos(z);

        Real df3_dx = -2 * std::sin(x) * std::cos(y) * std::cos(z);
        Real df3_dy = -2 * std::cos(x) * std::sin(y) * std::cos(z);
        Real df3_dz = -2 * std::cos(x) * std::cos(y) * std::sin(z);

        return f1 * df3_dx + f2 * df3_dy + f3 * df3_dz - (3 / Re) * f3;
    }
    
};




template<class OP>
void test_operator(VectorialFunction &function)
{
    const std::vector<size_t> N_values = {8, 16, 32, 64, 128};

    // Reynolds number
    const Real Re = 4700.0;

    // Define physical size of the problem (just for simplicity)
    const Real phy_dim = 1.0;


    for (const auto &N : N_values)
    {
        // Define dim as side dimension of the grid (just for simplicity)
        const index_t dim = N;

        // Define number of nodes for each axis
        const index_t nx = dim;
        const index_t ny = dim;
        const index_t nz = dim;

        // Define physical size of the problem for each axis
        const Real sx = phy_dim / (nx - 1);
        const Real sy = phy_dim / (ny - 1);
        const Real sz = phy_dim / (nz - 1);


        // Define initial velocity function
        auto initialVel = [](Real x, Real y, Real z) -> Vector {
            return {
                std::sin(x) * std::cos(y) * std::sin(z),
                std::cos(x) * std::sin(y) * std::sin(z),
                2 * std::cos(x) * std::cos(y) * std::cos(z)
            };
        };

        // Define initial pressure function
        // For the moment it does not work so do not care about it
        auto initialPres = [](Real x, Real y, Real z) -> Real {
            return 2.0;
        };
        

        // Define nodes
        std::array<index_t, 3> nodes = {nx, ny, nz};

        // Define spacing
        Vector spacing = {sx, sy, sz};

        // define the mesh:
        // instantiate from ghosted stagg grid
        // hint: second method, num of nodes in each direction, number of ghosts=1
        Grid<STANDARD> sg(nodes, spacing, 1);
        Grid<STANDARD> sg2(nodes, spacing, 1);

        sg.initGrid(initialVel, initialPres);
        set_boundary_values(sg, initialVel);


        // Test operators
        OP op(sg);
        IdentityOperator<STANDARD> I(sg);

        op.apply(sg2, sg);
        I.apply(sg, sg2);


        Real error = computeError(sg, function);
        std::cout << "h = " << phy_dim / dim << "\tN = " << dim << std::endl;
        std::cout << "\tError = " << error << std::endl;

    }
}





int main() {
    // try the solver

    // Reynolds number
    const Real Re = 4700.0;

    // Define dim as side dimension of the grid (just for simplicity)
    const index_t dim = 64;

    // Define number of nodes for each axis
    const index_t nx = dim;
    const index_t ny = dim;
    const index_t nz = dim;

    // Define physical size of the problem (just for simplicity)
    const Real phy_dim = 1.0;

    // Define physical size of the problem for each axis
    const Real sx = phy_dim / (nx - 1);
    const Real sy = phy_dim / (ny - 1);
    const Real sz = phy_dim / (nz - 1);

    // Define initial velocity function
    auto initialVel = [](Real x, Real y, Real z) -> Vector {
        return {
            std::sin(x) * std::cos(y) * std::sin(z),
            std::cos(x) * std::sin(y) * std::sin(z),
            2 * std::cos(x) * std::cos(y) * std::cos(z)
        };
    };

    // Define initial pressure function
    // For the moment it does not work so do not care about it
    auto initialPres = [](Real x, Real y, Real z) -> Real {
        return 2.0;
    };

    // Define nodes
    std::array<index_t, 3> nodes = {nx, ny, nz};

    // Define spacing
    Vector spacing = {sx, sy, sz};

    // define the mesh:
    // instantiate from ghosted stagg grid
    // hint: second method, num of nodes in each direction, number of ghosts=1
    Grid<STANDARD> sg(nodes, spacing, 1);
    Grid<STANDARD> sg2(nodes, spacing, 1);
    std::cout << "Grid created" << std::endl;

    // build the model
    // initialize the grid with initial values
    //Model<STANDARD> model(spacing, sg, Re, initialVel, initialPres);
    sg.initGrid(initialVel, initialPres);
    set_boundary_values(sg, initialVel);
    std::cout << "Model created, grid initialized" << std::endl;


    // ******************************************TESTS************************************************

    // To test an operator we can simply pass it as a template parameter to test_operator as follows
    // We have first to declare the real function to compare the operators
    
    // For the laplacian:
    {
        TestCaseLaplacian fun;
        std::cout << std::endl << "Testing the Laplacian Operator:" << std::endl;
        test_operator<LaplacianOperator<STANDARD>>(fun);
    }

    // For the convective term:
    {
        TestCaseConvection fun;
        std::cout << std::endl << "Testing the Convective Operator:" << std::endl;
        test_operator<ConvectiveOperator<STANDARD>>(fun);
    }


    // Alternatively it is possible to combine two or more operators in one, but I still have
    // to figure out how to generalize the testing functions

    LaplacianOperator<STANDARD> L(sg);
    IdentityOperator<STANDARD> I(sg);
    ConvectiveOperator<STANDARD> Unabla(sg);

    OneStepOperator<STANDARD> op = Unabla + (1 / Re) * L;

    // After defining the operator, just apply it on the grid
    op.apply(sg2, sg);

    // For now in order to measure the error we have to copy back the data from a grid to the other
    // TODO: try to figure out how to apply directly the operator to the grid
    I.apply(sg, sg2);

    TestCaseComposed fun(Re);
    //TestCaseIdentity fun;

    //set_values(model, fun);
    Real error = computeError(sg, fun);
    std::cout << std::endl << "Achieved error L2 norm for the composed operator:" << std::endl;
    std::cout << "Error = " << error << std::endl;
    std::cout << "N = " << dim << std::endl;



    // Output of last iteration
    /*
    VTKFile file = VTKConverter::exportModel(model, "testsolver output of last time iteration");
    file.writeFile("testsolver.vtk");
    */


    return 0;
}