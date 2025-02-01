#ifndef AEROHPC_A_PRINT_BUFFER_H
#define AEROHPC_A_PRINT_BUFFER_H

#include <fstream>
#include <string>
#include <iomanip>
#include <filesystem>
#include "Traits.hpp"

using namespace std;

using namespace filesystem;

class BufferPrinter {

    std::ofstream outputFile;
    const std::string printer_dir_path;

public:

    static constexpr int PRINT_VELOCITY = 1;
    static constexpr int PRINT_PRESSURE = 2;

    const int precision = 4;

    explicit BufferPrinter(const std::string &printer_dir_path) : printer_dir_path(printer_dir_path) {}

    void initDir() const {
#if DEBUG_PRINT_BUFFERS
        if (IS_MAIN_PROC) create_directories(printer_dir_path);
#endif

    }

    void setFile(const std::string &iterationFilePath) {
#if DEBUG_PRINT_BUFFERS
        if (outputFile.is_open()) outputFile.close();
        outputFile.open(printer_dir_path + "/" + iterationFilePath);
        outputFile << "y     z" << std::endl;
        outputFile << "↑   ↗" << std::endl;
        outputFile << "| ∕" << std::endl;
        outputFile << "0 — → x" << std::endl;
        outputFile << std::endl
        << "=============================================================================================="
        << "=============================================================================================="
        << "==============================================================================================" << std::endl;
#endif
    }

    void print(const Real *data, const std::string &buffer_name, const int type) {
        if (!outputFile.is_open()) return;

        outputFile << std::setprecision(precision) << std::fixed;

        outputFile << buffer_name << ":" << endl;
        std::string prev_space;
        std::string space;

        for (int k = -1; k <= params.loc_nZ ; ++k) prev_space += '\t';

        if (type & PRINT_VELOCITY) {
            for (int j = params.loc_nY; j >= -1; --j) {
                for (int k = params.loc_nZ; k >= -1; --k) {
                    outputFile << prev_space;
                    for (int i = -1; i <= params.loc_nX; ++i) {
                        outputFile << U(data, i, j, k) << " ";
                    }
                    outputFile << space << "\t" << prev_space;
                    for (int i = -1; i <= params.loc_nX; ++i) {
                        outputFile << V(data, i, j, k) << " ";
                    }
                    outputFile << space << "\t" << prev_space;
                    for (int i = -1; i <= params.loc_nX; ++i) {
                        outputFile << W(data, i, j, k) << " ";
                    }
                    outputFile << endl;
                    space += "\t";
                    prev_space = prev_space.substr(0, prev_space.size() - 1);
                }
                space = "";
                for (int k = -1; k <= params.loc_nZ ; ++k) prev_space += '\t';
                outputFile << endl;
            }
        }
        if (type & PRINT_PRESSURE) {
            for (int j = params.loc_nY; j >= -1; --j) {
                for (int k = -1; k <= params.loc_nZ ; ++k) {
                    outputFile << prev_space;
                    for (int i = -1; i <= params.loc_nX; ++i) {
                        outputFile << P(data, i, j, k) << " ";
                    }
                    outputFile << endl;
                    prev_space = prev_space.substr(0, prev_space.size() - 1);
                }
                space = "";
                for (int k = -1; k <= params.loc_nZ ; ++k) prev_space += '\t';
                outputFile << endl;
            }
        }
        outputFile << std::endl
        << "=============================================================================================="
        << "=============================================================================================="
        << "==============================================================================================" << std::endl;
    }
};

inline BufferPrinter bufferPrinter("./iterations");

#define enabledBufferPrinter if(IS_MAIN_PROC) bufferPrinter

#endif //AEROHPC_A_PRINT_BUFFER_H