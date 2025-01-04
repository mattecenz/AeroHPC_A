#ifndef AEROHPC_A_PRINT_BUFFER_H
#define AEROHPC_A_PRINT_BUFFER_H

#include "GridData.hpp"
#include <fstream>
#include <string>
#include <iomanip>
#include <filesystem>

using namespace std;

using namespace filesystem;

void print(GridData &grid, std::string &filename) {
    std::ofstream file;
    file.open(filename);

    file << std::setprecision(4) << std::fixed;

    int nx = grid.structure.nx;
    int ny = grid.structure.ny;
    int nz = grid.structure.nz;
    int gp = grid.structure.gp;

    file << "U:" << endl;
    std::string space;
    for (int j = ny - 1 + gp; j >= 0 - gp; --j) {
        for (int k = 0 - gp; k < nz + gp; ++k) {
            file << space;
            for (int i = 0 - gp; i < nx + gp; ++i) {
                file << grid.U(i, j, k) << " ";
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
                file << grid.V(i, j, k) << " ";
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
                file << grid.W(i, j, k) << " ";
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
                file << grid.P(i, j, k) << " ";
            }
            file << endl;
            space += "\t";
        }
        space = "";
        file << endl;
    }
}

#endif