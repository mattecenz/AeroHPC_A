#include "VTKConverter.hpp"

namespace VTKConverter {

    template<>
    std::vector<DataSection *> exportGrid(StaggeredGrid<STANDARD> &grid) {

        std::vector<std::vector<Real>> velocity;
        std::vector<std::vector<Real>> pressure;

        for (index_t z = 0; z < grid.nz; ++z) {
            for (index_t y = 0; y < grid.ny; ++y) {
                for (index_t x = 0; x < grid.nx; ++x) {
                    velocity.emplace_back(
                            std::vector<Real>{
                                    grid(U, x, y, z),
                                    grid(V, x, y, z),
                                    grid(W, x, y, z)
                            }
                    );
                    pressure.emplace_back(
                            std::vector<Real>{
                                    grid(P, x, y, z)
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

    template<>
    VTKFile exportModel(Model<STANDARD> &model, std::string description) {
        StaggeredGrid<STANDARD> grid = model.grid;

        VTKFile file({
                             static_cast<Real>(grid.nx) * model.dx,
                             static_cast<Real>(grid.ny) * model.dy,
                             static_cast<Real>(grid.nz) * model.dz,
                     }, {
                             static_cast<unsigned long>(grid.nx),
                             static_cast<unsigned long>(grid.ny),
                             static_cast<unsigned long>(grid.nz)
                     }, std::move(description));
        file << exportGrid(model.grid);

        return file;
    }
}