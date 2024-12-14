//
// Created by enti on 12/12/24.
//

#ifndef AEROHPC_A_JUSTFORTEST_HPP
#define AEROHPC_A_JUSTFORTEST_HPP


#include <iostream>
#include <fstream>
#include "GridData.hpp"

using namespace std;

void print(GridData &grid, std::string &filename) {

    std::ofstream file;
    file.open(filename);

    file << std::setprecision(2) << std::scientific;

    file << "U:," << endl;
    std::string space;
    for (int j = grid.structure.ny; j >= -1; --j) {
        for (int k = -1; k <= grid.structure.nz; ++k) {
            file << space;
            for (int i = -1; i <= grid.structure.nx; ++i) {
                file << grid.U(i, j, k) << ", ";
            }
            file << endl;
            space += "\t";
        }
        space = "";
        file << endl;
    }

    file << endl << "V:," << endl;
    for (int j = grid.structure.ny; j >= -1; --j) {
        for (int k = -1; k <= grid.structure.nz; ++k) {
            file << space;
            for (int i = -1; i <= grid.structure.nx; ++i) {
                file << grid.V(i, j, k) << ", ";
            }
            file << endl;
            space += "\t";
        }
        space = "";
        file << endl;
    }

    file << endl << "W:," << endl;
    for (int j = grid.structure.ny; j >= -1; --j) {
        for (int k = -1; k <= grid.structure.nz; ++k) {
            file << space;
            for (int i = -1; i <= grid.structure.nx; ++i) {
                file << grid.W(i, j, k) << ", ";
            }
            file << endl;
            space += "\t";
        }
        space = "";
        file << endl;
    }
}

void printBuffer(GridData &grid, std::string &filename) {

    std::ofstream file;
    file.open(filename);

    file << std::setprecision(2) << std::scientific;

    file << "U:" << endl;

    std::string space;
    for (int j = grid.structure.ny - 1; j > -1; --j) {
        for (int k = 0; k < grid.structure.nz; ++k) {
            file << space;
            for (int i = 0; i < grid.structure.nx; ++i) {
                file << grid.U(i, j, k) << " ";
            }
            file << endl;
            space += "\t";
        }
        space = "";
        file << endl;
    }

    file << endl << "V:" << endl;

    for (int j = grid.structure.ny - 1; j > -1; --j) {
        for (int k = 0; k < grid.structure.nz; ++k) {
            file << space;
            for (int i = 0; i < grid.structure.nx; ++i) {
                file << grid.V(i, j, k) << " ";
            }
            file << endl;
            space += "\t";
        }
        space = "";
        file << endl;
    }

    file << endl << "W:" << endl;

    for (int j = grid.structure.ny - 1; j > -1; --j) {
        for (int k = 0; k < grid.structure.nz; ++k) {
            file << space;
            for (int i = 0; i < grid.structure.nx; ++i) {
                file << grid.W(i, j, k) << " ";
            }
            file << endl;
            space += "\t";
        }
        space = "";
        file << endl;
    }

}



#endif //AEROHPC_A_JUSTFORTEST_HPP
