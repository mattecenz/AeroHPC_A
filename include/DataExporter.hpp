#ifndef DATAEXPORTER_HPP
#define DATAEXPORTER_HPP

#include "Traits.hpp"
#include "data/SolverData.hpp"
#include <mpi.h>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <iomanip>
#include "vtk/VTKFile.hpp"

inline void probePoint(std::vector<Real> &velocity, std::vector<Real> &pressure, const Real x, const Real y, const Real z) {
    const Real local_x = x - params.originX;
    const Real local_y = y - params.originY;
    const Real local_z = z - params.originZ;

    const auto prev_x = static_cast<index_t>(floor(local_x / params.dX)) - params.st_nX;
    const auto next_x = static_cast<index_t>(ceil(local_x / params.dX)) - params.st_nX;

    const auto prev_y = static_cast<index_t>(floor(local_y / params.dY)) - params.st_nY;
    const auto next_y = static_cast<index_t>(ceil(local_y / params.dY)) - params.st_nY;

    const auto prev_z = static_cast<index_t>(floor(local_z / params.dZ)) - params.st_nZ;
    const auto next_z = static_cast<index_t>(ceil(local_z / params.dZ)) - params.st_nZ;

    const Real prev_x_perc = 1.0 - ((local_x - prev_x) / params.dX);
    const Real next_x_perc = 1.0 - prev_x_perc;

    const Real prev_y_perc = 1.0 - ((local_y - prev_y) / params.dY);
    const Real next_y_perc = 1.0 - prev_y_perc;

    const Real prev_z_perc = 1.0 - ((local_z - prev_z) / params.dZ);
    const Real next_z_perc = 1.0 - prev_z_perc;

#define getComponentLetter(letter, x, y, z) \
    const Real letter##U = PU(interpData_ptr, x, y, z); \
    const Real letter##V = PV(interpData_ptr, x, y, z); \
    const Real letter##W = PW(interpData_ptr, x, y, z); \
    const Real letter##P = PP(interpData_ptr, x, y, z)

    getComponentLetter(A, prev_x, prev_y, next_z);
    getComponentLetter(B, next_x, prev_y, next_z);
    getComponentLetter(C, next_x, next_y, next_z);
    getComponentLetter(D, prev_x, next_y, next_z);

    getComponentLetter(E, prev_x, next_y, prev_z);
    getComponentLetter(F, prev_x, prev_y, prev_z);
    getComponentLetter(G, next_x, prev_y, prev_z);
    getComponentLetter(H, next_x, next_y, prev_z);
#undef getComponentLetter

    const Real px = 1 + prev_x_perc;
    const Real nx = 1 + next_x_perc;
    const Real py = 1 + prev_y_perc;
    const Real ny = 1 + next_y_perc;
    const Real pz = 1 + prev_z_perc;
    const Real nz = 1 + next_z_perc;

#define getInterpolatedLetter(letter, perc_prev, perc_next, prev_letter, next_letter) \
    const Real letter##U = perc_prev * prev_letter##U + perc_next * next_letter##U; \
    const Real letter##V = perc_prev * prev_letter##V + perc_next * next_letter##V; \
    const Real letter##W = perc_prev * prev_letter##W + perc_next * next_letter##W; \
    const Real letter##P = perc_prev * prev_letter##P + perc_next * next_letter##P

    // weighted interpolation along y dir
    getInterpolatedLetter(J, py, ny, A, D);
    getInterpolatedLetter(K, py, ny, B, C);
    getInterpolatedLetter(L, py, ny, G, H);
    getInterpolatedLetter(M, py, ny, F, E);

    // weighted interpolation along x dir
    getInterpolatedLetter(U, px, nx, J, K);
    getInterpolatedLetter(O, px, nx, M, L);

    // weighted interpolation along z dir
    getInterpolatedLetter(probe, pz, nz, O, U);

#undef getInterpolatedLetter

    // add normalized values
    velocity.push_back(probeU / 27.0);
    velocity.push_back(probeV / 27.0);
    velocity.push_back(probeW / 27.0);
    pressure.push_back(probeP / 27.0);
}

inline void extractFaceData(std::vector<Real> &point_coord,
                            std::vector<Real> &velocity_data,
                            std::vector<Real> &pressure_data,
                            const std::array<Real, 3> &point) {
    const Real dx = params.dX;
    const Real dy = params.dY;
    const Real dz = params.dZ;

    const index_t nx = params.phy_nX - 1;
    const index_t ny = params.phy_nY - 1;
    const index_t nz = params.phy_nZ - 1;

    /// Get point index if it is in this domain
    const index_t px = params.st_nX;
    const index_t py = params.st_nY;
    const index_t pz = params.st_nZ;

    const Real origin_x = real(px) * dx + params.originX;
    const Real origin_y = real(py) * dy + params.originY + (!params.isOnBottom * dy);
    const Real origin_z = real(pz) * dz + params.originZ + (!params.isOnLeft * dz);

    const Real end_x = origin_x + real(nx) * dx;
    const Real end_y = origin_y + real(ny) * dy;
    const Real end_z = origin_z + real(nz) * dz;

    const bool is_in_x_range = (point[0] >= origin_x && point[0] <= end_x);
    const bool is_in_y_range = (point[1] >= origin_y && point[1] <= end_y);
    const bool is_in_z_range = (point[2] >= origin_z && point[2] <= end_z);


    if (is_in_x_range) {
        for (Real y = origin_y; y <= end_y; y += dy) {
            for (Real z = origin_z; z <= end_z; z += dz) {
                point_coord.push_back(point[0]);
                point_coord.push_back(y);
                point_coord.push_back(z);
                probePoint(velocity_data, pressure_data, point[0], y, z);
            }
        }
    }

    if (is_in_y_range) {
        for (Real x = origin_x; x <= end_x; x += dx) {
            for (Real z = origin_z; z <= end_z; z += dz) {
                point_coord.push_back(x);
                point_coord.push_back(point[1]);
                point_coord.push_back(z);
                probePoint(velocity_data, pressure_data, x, point[1], z);
            }
        }
    }

    if (is_in_z_range) {
        for (Real x = origin_x; x <= end_x; x += dx) {
            for (Real y = origin_y; y <= end_y; y += dy) {
                point_coord.push_back(x);
                point_coord.push_back(y);
                point_coord.push_back(point[2]);
                probePoint(velocity_data, pressure_data, x, y, point[2]);
            }
        }
    }
}


inline void extractLineData(std::vector<Real> &point_coord,
                            std::vector<Real> &velocity_data,
                            std::vector<Real> &pressure_data,
                            const int axis,
                            const std::array<Real, 3> &point) {

    const Real dx = params.dX;
    const Real dy = params.dY;
    const Real dz = params.dZ;

    const index_t nx = params.phy_nX - 1;
    const index_t ny = params.phy_nY - 1;
    const index_t nz = params.phy_nZ - 1;

    /// Get point index if it is in this domain
    const index_t px = params.st_nX;
    const index_t py = params.st_nY;
    const index_t pz = params.st_nZ;

    const Real origin_x = real(px) * dx + params.originX;
    const Real origin_y = real(py) * dy + params.originY + (!params.isOnBottom * dy);
    const Real origin_z = real(pz) * dz + params.originZ + (!params.isOnLeft * dz);

    const Real end_x = origin_x + real(nx) * dx;
    const Real end_y = origin_y + real(ny) * dy;
    const Real end_z = origin_z + real(nz) * dz;

    const bool is_in_x_range = (point[0] >= origin_x && point[0] <= end_x);
    const bool is_in_y_range = (point[1] >= origin_y && point[1] <= end_y);
    const bool is_in_z_range = (point[2] >= origin_z && point[2] <= end_z);

    // If ask for x-axis then the point should be in range of y and z
    if (axis == 0 && is_in_y_range && is_in_z_range) {
        for (Real x = origin_x; x <= end_x; x += dx) {
            point_coord.push_back(x);
            point_coord.push_back(point[1]);
            point_coord.push_back(point[2]);
            probePoint(velocity_data, pressure_data, x, point[1], point[2]);
        }
    }

    // If ask for y-axis then the point should be in range of x and z
    if (axis == 1 && is_in_x_range && is_in_z_range) {
        for (Real y = origin_y; y <= end_y; y += dy) {
            point_coord.push_back(point[0]);
            point_coord.push_back(y);
            point_coord.push_back(point[2]);
            probePoint(velocity_data, pressure_data, point[0], y, point[2]);
        }
    }

    // If ask for z-axis then the point should be in range of x and y
    if (axis == 2 && is_in_x_range && is_in_y_range) {
        for (Real z = origin_z; z <= end_z; z += dz) {
            point_coord.push_back(point[0]);
            point_coord.push_back(point[1]);
            point_coord.push_back(z);
            probePoint(velocity_data, pressure_data, point[0], point[1], z);
        }
    }
}


inline void writeVtkFile(
    const std::string &filename,
    const std::string &description,
    std::vector<Real> &point_values,
    std::vector<Real> &velocity_values,
    std::vector<Real> &pressure_values) {
    const int local_n = static_cast<int>(point_values.size() / 3);

    for (Real &point_value: point_values)
        swap_endian(point_value);

    for (Real &velocity_value: velocity_values)
        swap_endian(velocity_value);

    for (Real &pressure_value: pressure_values)
        swap_endian(pressure_value);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_File file;
    MPI_Status status;

    constexpr unsigned long datasize = sizeof(Real);
    const MPI_Datatype mpi_datatype = Real_MPI;
    const std::string &type = Real_Dataname;

    int n = 0;
    MPI_Reduce(&local_n, &n, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast n to all processes

    // Calculate ASCII header and section sizes (constant across all processes)
    std::string header = "# vtk DataFile Version 2.0\n" + description + "\nBINARY\nDATASET UNSTRUCTURED_GRID\n";
    std::string points_header = "POINTS " + std::to_string(n) + " " + type + "\n";
    std::string velocity_header = "\nPOINT_DATA " + std::to_string(n) + "\nSCALARS U " + type + " 3\nLOOKUP_TABLE default\n";
    std::string pressure_header = "\nSCALARS P " + type + " 1\nLOOKUP_TABLE default\n";

    MPI_Offset offset = 0;
    MPI_Offset temp_local_offset = 0;

    // Open the file for writing
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

    // WRITE VTK HEADER AND POINTS HEADER ////////////////////////////////////////////////////////////////////////////////
    if (rank == 0) {
        // Step 1: Write ASCII header
        MPI_File_write_at(file, 0, header.c_str(), header.size(), MPI_CHAR, &status);

        // Step 2: Write ASCII points section
        MPI_File_write_at(file, header.size(), points_header.c_str(), points_header.size(), MPI_CHAR, &status);
    }

    // update offset (vtk_h + pnt_h)
    offset += header.size() + points_header.size();

    MPI_Barrier(MPI_COMM_WORLD);

    // WRITE POINT SECTION //////////////////////////////////////////////////////////////////////////////////////////////
    temp_local_offset = 0;
    MPI_Exscan(&local_n, &temp_local_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) temp_local_offset = 0;

    temp_local_offset = temp_local_offset * datasize * 3;

    MPI_File_write_at(file, offset + temp_local_offset, point_values.data(), local_n * 3, mpi_datatype, &status);

    // update offset (vtk_h + pnt_h + pnt_sec)
    offset += n * 3 * datasize;

    MPI_Barrier(MPI_COMM_WORLD);

    // WRITE VELOCITY HEADER ////////////////////////////////////////////////////////////////////////////////////////////
    if (rank == 0) {
        // Step 4: Write ASCII scalars section
        MPI_File_write_at(file, offset, velocity_header.c_str(), velocity_header.size(), MPI_CHAR, &status);
    }

    // update offset (vtk_h + pnt_h + pnt_sec + vel_h)
    offset += velocity_header.size();

    MPI_Barrier(MPI_COMM_WORLD);

    // WRITE VELOCITY SECTION ///////////////////////////////////////////////////////////////////////////////////////////
    temp_local_offset = 0;
    MPI_Exscan(&local_n, &temp_local_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) temp_local_offset = 0;

    temp_local_offset = temp_local_offset * datasize * 3;

    MPI_File_write_at(file, offset + temp_local_offset, velocity_values.data(), local_n * 3, mpi_datatype, &status);

    // update offset (vtk_h + pnt_h + pnt_sec + vel_h + vel_sec)
    offset += n * 3 * datasize;

    MPI_Barrier(MPI_COMM_WORLD);

    // WRITE PRESSURE HEADER ////////////////////////////////////////////////////////////////////////////////////////////
    if (rank == 0) {
        // Step 4: Write ASCII scalars section
        MPI_File_write_at(file, offset, pressure_header.c_str(), pressure_header.size(), MPI_CHAR, &status);
    }

    // update offset (vtk_h + pnt_h + pnt_sec + vel_h + vel_sec + pres_h)
    offset += pressure_header.size();

    MPI_Barrier(MPI_COMM_WORLD);

    // WRITE PRESSURE SECTION ///////////////////////////////////////////////////////////////////////////////////////////
    temp_local_offset = 0;
    MPI_Exscan(&local_n, &temp_local_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) temp_local_offset = 0;

    temp_local_offset = temp_local_offset * datasize;

    MPI_File_write_at(file, offset + temp_local_offset, pressure_values.data(), local_n, mpi_datatype, &status);

    MPI_Barrier(MPI_COMM_WORLD);

    // Close the file
    MPI_File_close(&file);
}

inline void writeDatFile(
    const std::string &filename,
    const std::vector<Real> &point_values,
    const std::vector<Real> &velocity_values,
    const std::vector<Real> &pressure_values) {
    const int local_n = static_cast<int>(point_values.size() / 3);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_File file;
    MPI_Status status;

    int n = 0;
    MPI_Reduce(&local_n, &n, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast n to all processes

    // Calculate ASCII header and section sizes (constant across all processes)
    MPI_Offset offset = 0;

    // Open the file for writing
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

    std::stringstream ss;
    ss << std::scientific << std::setprecision(8);

    for (int i = 0; i < local_n; i++) {
        ss << point_values[i * 3] << " " << point_values[i * 3 + 1] << " " << point_values[i * 3 + 2];
        ss << " " << velocity_values[i * 3] << " " << velocity_values[i * 3 + 1] << " " << velocity_values[i * 3 + 2];
        ss << " " << pressure_values[i * 3] << "\n";
    }

    const std::string str = ss.str();
    const size_t size = str.size();
    // WRITE POINT SECTION //////////////////////////////////////////////////////////////////////////////////////////////
    MPI_Exscan(&size, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) offset = 0;

    MPI_File_write_at(file, offset, str.c_str(), str.size(), MPI_CHAR, &status);
    MPI_Barrier(MPI_COMM_WORLD);

    // Close the file
    MPI_File_close(&file);
}

#endif //DATAEXPORTER_HPP
