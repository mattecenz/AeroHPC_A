#include "vtk/VTKConverter.hpp"

namespace VTKConverter
{

    std::vector<DataSection *> exportData(Real *grid)
    {

        std::vector<std::vector<Real>> velocity;
        std::vector<std::vector<Real>> pressure;

        for (index_t z = -1; z < params.loc_nZ+1; ++z)
        {
            for (index_t y = -1; y < params.loc_nY+1; ++y)
            {
                for (index_t x = -1; x < params.loc_nX+1; ++x)
                {
                    velocity.emplace_back(
                        std::vector<Real>{
                            U(grid, x, y, z),
                            V(grid, x, y, z),
                            W(grid, x, y, z)});
                    pressure.emplace_back(
                        std::vector<Real>{
                            P(grid, x, y, z)});
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

    VTKFile exportGrid(Real *grid, std::string description)
    {
        VTKFile file({real(params.loc_nX+2) * params.dX,
                      real(params.loc_nY+2) * params.dY,
                      real(params.loc_nZ+2) * params.dZ},
                     {static_cast<unsigned long>(params.loc_nX+2),
                      static_cast<unsigned long>(params.loc_nY+2),
                      static_cast<unsigned long>(params.loc_nZ+2)},
                     std::move(description));
        file << exportData(grid);

        return file;
    }
}