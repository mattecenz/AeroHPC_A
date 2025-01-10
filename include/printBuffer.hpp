#ifndef AEROHPC_A_PRINT_BUFFER_H
#define AEROHPC_A_PRINT_BUFFER_H

#include "GridData.hpp"
#include <fstream>
#include <string>
#include <iomanip>
#include <filesystem>

using namespace std;

using namespace filesystem;

void print(GridData &grid, std::string &filename, const std::string &name, bool print_vel, bool print_P) {
    std::ofstream file;
    file.open(filename,ios_base::out | ios_base::app);

    file << std::setprecision(6) << std::fixed;

    int nx = grid.structure.nx;
    int ny = grid.structure.ny;
    int nz = grid.structure.nz;
    int gp = grid.structure.gp;

    std::string space;

    if (print_vel) {
        file << name << ":" << endl;
        for (int j = ny - 1 + gp; j >= 0 - gp; --j) {
            for (int k = 0 - gp; k < nz + gp; ++k) {
                file << space;
                for (int i = 0 - gp; i < nx + gp; ++i) {
                    file << "( " << grid.U(i, j, k) << " , " << grid.V(i,j,k) << ", " << grid.W(i,j,k) << " )" << " ";
                }
                file << endl;
                space += "\t";
            }
            space = "";
            file << endl;
        }
    }
    if (print_P) {
        file << endl << name <<":" << endl;
        for (int j = ny - 1 + gp; j >= 0 - gp; --j) {
            for (int k = 0 - gp; k < nz + gp; ++k) {
                file << space;
                for (int i = 0 - gp; i < nx + gp; ++i) {
                    file << grid.P(i, j, k) << " ";
                }
                file << endl;
                space += "\t";
            }
            space = "";
            file << endl;
        }
    }
}

#ifdef DEBUG_PRINT_BUFFERS
#define b_print(buff, dir, n, name, p_Vel, p_P) \
if (!rank) { \
std::string nn{dir + to_string(n)}; \
print(buff, nn, name, p_Vel, p_P); \
}
#define c_dir(dir) \
if (!rank) {\
create_directories(dir); \
}
#else
#define b_print(a,b,c,d,p_Vel, p_P) //
#define c_dir(a) //
#endif

#endif