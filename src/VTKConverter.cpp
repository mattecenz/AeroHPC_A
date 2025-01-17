#include "vtk/VTKConverter.hpp"
#include "data/SolverData.hpp"

namespace VTKConverter
{

    std::vector<DataSection *> exportData(const Real *physical_grid)
    {
        std::vector<std::vector<Real>> velocity;
        std::vector<std::vector<Real>> pressure;

        for (index_t z = 0; z < params.phy_nZ; ++z)
        {
            for (index_t y = 0; y < params.phy_nY; ++y)
            {
                for (index_t x = 0; x < params.phy_nX; ++x)
                {
                    velocity.emplace_back(
                        std::vector<Real>{
                            PU(physical_grid, x, y, z),
                            PV(physical_grid, x, y, z),
                            PW(physical_grid, x, y, z)});
                    pressure.emplace_back(
                        std::vector<Real>{
                            PP(physical_grid, x, y, z)});
                }
            }
        }
        DataSection *velocity_section = new ScalarsSection<Real>("velocity", velocity);
        DataSection *pressure_section = new ScalarsSection<Real>("pressure", pressure);
        std::vector<DataSection *> out;
        out.push_back(velocity_section);
        out.push_back(pressure_section);
        return out;
    }

    VTKFile exportGrid(const Real *grid, std::string description)
    {
        VTKFile file({params.dX, params.dY, params.dZ},
                   {params.originX, params.originY, params.originZ},
                     {static_cast<unsigned long>(params.phy_nX),
                      static_cast<unsigned long>(params.phy_nY),
                      static_cast<unsigned long>(params.phy_nZ)},
                     std::move(description));
        file << exportData(grid);

        return file;
    }
}