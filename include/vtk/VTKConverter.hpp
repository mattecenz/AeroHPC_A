#ifndef AEROHPC_A_VTKCONVERTER_H
#define AEROHPC_A_VTKCONVERTER_H

#include "vtk/VTKFile.hpp"
#include "data/SolverData.hpp"

namespace VTKConverter{

    VTKFile exportGrid(const Real* grid, std::string description = "");

    std::vector<DataSection *> exportData(const Real* physical_grid);
}

#endif //AEROHPC_A_VTKCONVERTER_H
