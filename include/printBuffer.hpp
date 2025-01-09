#ifndef AEROHPC_A_PRINT_BUFFER_H
#define AEROHPC_A_PRINT_BUFFER_H

#include <fstream>
#include <string>
#include <iomanip>
#include <filesystem>
#include "SolverData.hpp"

using namespace std;

using namespace filesystem;

inline void print(const Real *data, std::string &filename, bool has_ghosts) {
    std::ofstream file;
    file.open(filename);

    file << std::setprecision(4) << std::fixed;

    int nx = has_ghosts ? params.loc_gnX : params.loc_nX;
    int ny = has_ghosts ? params.loc_gnY : params.loc_nY;
    int nz = has_ghosts ? params.loc_gnZ : params.loc_nZ;
    int gp = has_ghosts ? 1 : 0;

    file << "U:" << endl;
    std::string space;
    for (int j = ny - 1 + gp; j >= 0 - gp; --j) {
        for (int k = 0 - gp; k < nz + gp; ++k) {
            file << space;
            for (int i = 0 - gp; i < nx + gp; ++i) {
                file << U(data, i, j, k) << " ";
            }
            file << endl;
            space += "\t";
        }
        space = "";
        file << endl;
    }

    file << endl << "V:" << endl;
    for (int j = ny - 1 + gp; j >= 0 - gp; --j) {
        for (int k = 0 - gp; k < nz + gp; ++k) {
            file << space;
            for (int i = 0 - gp; i < nx + gp; ++i) {
                file << V(data, i, j, k) << " ";
            }
            file << endl;
            space += "\t";
        }
        space = "";
        file << endl;
    }

    file << endl << "W:" << endl;
    for (int j = ny - 1 + gp; j >= 0 - gp; --j) {
        for (int k = 0 - gp; k < nz + gp; ++k) {
            file << space;
            for (int i = 0 - gp; i < nx + gp; ++i) {
                file << W(data, i, j, k) << " ";
            }
            file << endl;
            space += "\t";
        }
        space = "";
        file << endl;
    }

    file << endl << "P:" << endl;
    for (int j = ny - 1 + gp; j >= 0 - gp; --j) {
        for (int k = 0 - gp; k < nz + gp; ++k) {
            file << space;
            for (int i = 0 - gp; i < nx + gp; ++i) {
                file << P(data, i, j, k) << " ";
            }
            file << endl;
            space += "\t";
        }
        space = "";
        file << endl;
    }
}

inline std::string dir;

#ifdef DEBUG_PRINT_BUFFERS
#define b_print(buff, n, has_ghost) \
if (!rank) { \
std::string nn{dir + to_string(n)}; \
print(buff, nn, has_ghost); \
}
#define c_dir() \
if (!rank) {\
create_directories(dir); \
}
#else
#define b_print(a,b,c) //
#define c_dir(a) //
#endif

#endif