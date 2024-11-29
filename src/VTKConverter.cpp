#include "VTKConverter.hpp"

namespace VTKConverter {

    std::vector<DataSection *> exportData(GridData &grid) {

        std::vector<std::vector<Real>> velocity;
        std::vector<std::vector<Real>> pressure;

        for (index_t z = 0; z < grid.structure.nz; ++z) {
            for (index_t y = 0; y < grid.structure.ny; ++y) {
                for (index_t x = 0; x < grid.structure.nx; ++x) {
                    velocity.emplace_back(
                            std::vector<Real>{
                                    grid.U(x, y, z),
                                    grid.V(x, y, z),
                                    grid.W(x, y, z)
                            }
                    );
                    pressure.emplace_back(
                            std::vector<Real>{
                                    grid.P(x, y, z)
                            }
                    );
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

    VTKFile exportGrid(GridData &grid, std::string description) {
        VTKFile file({
                             real(grid.structure.nx) * grid.structure.dx,
                             real(grid.structure.ny) * grid.structure.dy,
                             real(grid.structure.nz) * grid.structure.dz,
                     }, {
                             static_cast<unsigned long>(grid.structure.nx),
                             static_cast<unsigned long>(grid.structure.ny),
                             static_cast<unsigned long>(grid.structure.nz)
                     }, std::move(description));
        file << exportData(grid);

        return file;
    }
}