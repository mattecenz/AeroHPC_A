#ifndef AEROHPC_A_VTKCONVERTER_H
#define AEROHPC_A_VTKCONVERTER_H

#include "VTKFile.hpp"
#include "Model.hpp"
#include "StaggeredGrid.hpp"

namespace VTKConverter{

    template <Addressing_T A>
    VTKFile exportModel(Model<A>& model, std::string description = "");

    template <Addressing_T A>
    std::vector<DataSection *> exportGrid(StaggeredGrid<A>& grid);
}

#endif //AEROHPC_A_VTKCONVERTER_H
