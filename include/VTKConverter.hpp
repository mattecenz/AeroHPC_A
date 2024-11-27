#ifndef AEROHPC_A_VTKCONVERTER_H
#define AEROHPC_A_VTKCONVERTER_H

#include "VTKFile.hpp"
#include "Grid.hpp"

namespace VTKConverter{

    VTKFile exportGrid(Grid& grid, std::string description = "");

    std::vector<DataSection *> exportData(Grid& grid);
}

#endif //AEROHPC_A_VTKCONVERTER_H
