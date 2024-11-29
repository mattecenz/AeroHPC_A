#ifndef AEROHPC_A_VTKCONVERTER_H
#define AEROHPC_A_VTKCONVERTER_H

#include "VTKFile.hpp"
#include "GridData.hpp"

namespace VTKConverter{

    VTKFile exportGrid(GridData& grid, std::string description = "");

    std::vector<DataSection *> exportData(GridData& grid);
}

#endif //AEROHPC_A_VTKCONVERTER_H
