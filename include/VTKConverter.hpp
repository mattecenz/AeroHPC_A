#ifndef AEROHPC_A_VTKCONVERTER_H
#define AEROHPC_A_VTKCONVERTER_H

#include "VTKFile.hpp"
#include "Grid.hpp"

namespace VTKConverter{

    template <Addressing_T A>
    VTKFile exportGrid(Grid<A>& grid, std::string description = "");

    template <Addressing_T A>
    std::vector<DataSection *> exportData(Grid<A>& grid);
}

#endif //AEROHPC_A_VTKCONVERTER_H
