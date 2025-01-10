#include "Traits.hpp"
#include "mathUtils.hpp"
#include "Logger.hpp"
#include "Boundaries.hpp"
#include <fstream>
#include <cmath>
#include <iomanip>
#include <iostream>

class TestFunction {
    public:
        static Real u(Real x, Real y, Real z){
            return std::sin(x)+std::sin(y)+std::sin(z);
        }
    
        static Real v(Real x, Real y, Real z){
            return std::cos(x)+std::cos(y)+std::cos(z);
        }
    
        static Real w(Real x, Real y, Real z){
            return x*x+y*y+z*z;
        }

        static Real div(Real x, Real y, Real z){
            return std::cos(x)-std::sin(y)+2*z;
        }
};    

int main(int argc, char **argv)
{
    bool writeToFile = false;
    std::ofstream outputFile("results.txt");
    const int npy = 1;
    const int npz = 1;

    const std::vector<int> nodes = {8,16,32,64,128};
    const Real dim_x = 1.0;
    const Real dim_y = 1.0;
    const Real dim_z = 1.0;

    logger.openSection("Running the Divergence test").spacer().printValue(5, "Phy dim", std::to_string(dim_x) + " x " + std::to_string(dim_y) + " x " + std::to_string(dim_z));

    auto initialVel = [](Real x, Real y, Real z) -> Vector {
        return {TestFunction::u(x, y, z), TestFunction::v(x, y, z), TestFunction::w(x, y, z)};
    };

    // Define initial pressure function (not used, but required)
    auto initialPres = [](Real x, Real y, Real z) -> Real {
        return 0.0;
    };

    

    std::vector<Real> error;

    logger.spacer().openTable("Number of nodes",{"l2error"});

    for (int n : nodes)
    {
        const Real sx = dim_x / real(n);
        const Real sy = dim_y / real(n);
        const Real sz = dim_z / real(n);
        Vector spacing = {sx, sy, sz};
        Idx3 nodes = {n, n, n};
        Idx3 displacement = {0, 0, 0};

        GridStructure modelStructure(nodes, spacing, displacement, 1);
        GridData input(modelStructure);

        // (Manually) Set all points, including ghosts (this are set automatically in the whole thing) (allegedly)
        for(int i=-1;i<n+1;i++){
            for(int j=-1;j<n+1;j++){
                for(int k=-1;k<n+1;k++){
                    Real x = real(i) * input.structure.dx;
                    Real y = real(j) * input.structure.dy;
                    Real z = real(k) * input.structure.dz;
                    input.U(i,j,k) = TestFunction::u(x,y,z);
                    input.V(i,j,k) = TestFunction::v(x,y,z);
                    input.W(i,j,k) = TestFunction::w(x,y,z);
                }
            }
        }


        //Compute divergence and error
        Real sum = 0.0;

        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                for(int k=0;k<n;k++){
                    // Convert grid indices to real space coordinates
                    Real x = real(i) * input.structure.dx;
                    Real y = real(j) * input.structure.dy;
                    Real z = real(k) * input.structure.dz;

                    // Calculate the difference
                    Real diff = mathUtils::vel_div(input,i,j,k) - TestFunction::div(x,y,z);

                    if(writeToFile) outputFile  << std::fixed << std::setprecision(6) << diff << "\t";

                    // Add the squares of the differences to sum
                    sum += (diff * diff);
                }
                if(writeToFile) outputFile << std::endl;
            }
            if(writeToFile) outputFile << std::endl;
        }

        Real l2norm = sum * input.structure.dx * input.structure.dy * input.structure.dz;
        error.push_back(l2norm);

        logger.printTableValues(n,{l2norm});
        outputFile << std::endl;
    }

    logger.closeTable();
    logger.closeSection();

    outputFile.close();

    //Write to csv for plotting script
    std::ofstream csvFile("output.csv");
    csvFile << "size,error" << std::endl;
    for (int i = 0; i < nodes.size(); ++i) csvFile << nodes[i] << "," << error[i] << std::endl;

    return 0;
}
