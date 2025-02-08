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
#include <iostream>

inline void probePoint(std::vector<Real> &velocity, std::vector<Real> &pressure, const Real x, const Real y, const Real z) {
    // center the domain onto cartesian grid
    const Real centered_x = x - params.originX;
    const Real centered_y = y - params.originY;
    const Real centered_z = z - params.originZ;

    // get coordinates on subdomain system
    const Real local_x = centered_x - (params.st_nX * params.dX);
    const Real local_y = centered_y - (params.st_nY * params.dY);
    const Real local_z = centered_z - (params.st_nZ * params.dZ);

    // get nodes on subdomain system
    const auto prev_i = static_cast<index_t>(floor(local_x / params.dX));
    const auto next_i = static_cast<index_t>(ceil(local_x / params.dX));

    const auto prev_j = static_cast<index_t>(floor(local_y / params.dY));
    const auto next_j = static_cast<index_t>(ceil(local_y / params.dY));

    const auto prev_k = static_cast<index_t>(floor(local_z / params.dZ));
    const auto next_k = static_cast<index_t>(ceil(local_z / params.dZ));

    if (prev_i < 0) std::cout << "[ERROR] prev_i = " << prev_i << std::endl;
    if (prev_j < 0) std::cout << "[ERROR] prev_j = " << prev_j << std::endl;
    if (prev_k < 0) std::cout << "[ERROR] prev_k = " << prev_k << std::endl;
    if (next_i > params.phy_nX - 1) std::cout << "[ERROR] next_i = " << next_i << std::endl;
    if (next_j > params.phy_nY - 1) std::cout << "[ERROR] next_j = " << next_j << std::endl;
    if (next_k > params.phy_nZ - 1) std::cout << "[ERROR] next_k = " << next_k << std::endl;

    // get coordinates of nodes
    const Real prev_x = real(prev_i) * params.dX;
    const Real prev_y = real(prev_j) * params.dY;
    const Real prev_z = real(prev_k) * params.dZ;

    // get ratio of contribution
    const Real pprev_x = 1.0 - ((local_x - prev_x) / params.dX);
    const Real pnext_x = 1.0 - pprev_x;

    const Real pprev_y = 1.0 - ((local_y - prev_y) / params.dY);
    const Real pnext_y = 1.0 - pprev_y;

    const Real pprev_z = 1.0 - ((local_z - prev_z) / params.dZ);
    const Real pnext_z = 1.0 - pprev_z;

#define getComponentLetter(letter, i, j, k) \
    const Real letter##U = PU(interpData_ptr, i, j, k); \
    const Real letter##V = PV(interpData_ptr, i, j, k); \
    const Real letter##W = PW(interpData_ptr, i, j, k); \
    const Real letter##P = PP(interpData_ptr, i, j, k)

    getComponentLetter(A, prev_i, prev_j, next_k);
    getComponentLetter(B, next_i, prev_j, next_k);
    getComponentLetter(C, next_i, next_j, next_k);
    getComponentLetter(D, prev_i, next_j, next_k);

    getComponentLetter(E, prev_i, next_j, prev_k);
    getComponentLetter(F, prev_i, prev_j, prev_k);
    getComponentLetter(G, next_i, prev_j, prev_k);
    getComponentLetter(H, next_i, next_j, prev_k);
#undef getComponentLetter

#define getInterpolatedLetter(letter, perc_prev, perc_next, prev_letter, next_letter) \
    const Real letter##U = perc_prev * prev_letter##U + perc_next * next_letter##U; \
    const Real letter##V = perc_prev * prev_letter##V + perc_next * next_letter##V; \
    const Real letter##W = perc_prev * prev_letter##W + perc_next * next_letter##W; \
    const Real letter##P = perc_prev * prev_letter##P + perc_next * next_letter##P

    // weighted interpolation along y dir
    getInterpolatedLetter(J, pprev_y, pnext_y, A, D);
    getInterpolatedLetter(K, pprev_y, pnext_y, B, C);
    getInterpolatedLetter(L, pprev_y, pnext_y, G, H);
    getInterpolatedLetter(M, pprev_y, pnext_y, F, E);

    // weighted interpolation along x dir
    getInterpolatedLetter(U, pprev_x, pnext_x, J, K);
    getInterpolatedLetter(O, pprev_x, pnext_x, M, L);

    // weighted interpolation along z dir
    getInterpolatedLetter(probe, pprev_z, pnext_z, O, U);

#undef getInterpolatedLetter

    // add normalized values
    velocity.push_back(probeU);
    velocity.push_back(probeV);
    velocity.push_back(probeW);
    pressure.push_back(probeP);
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

    //                   domain origin  +   subdomain start       + optional padding for subdomain intersection discrimination
    const Real start_x = params.originX + real(params.st_nX) * dx;
    const Real start_y = params.originY + real(params.st_nY) * dy + (!params.isOnBottom * dy);
    const Real start_z = params.originZ + real(params.st_nZ) * dz + (!params.isOnLeft * dz);

    //                  domain origin +     subdomain start     + subdomain size
    const Real end_x = params.originX + real(params.st_nX) * dx + real(nx) * dx;
    const Real end_y = params.originY + real(params.st_nY) * dy + real(ny) * dy;
    const Real end_z = params.originZ + real(params.st_nZ) * dz + real(nz) * dz;

    const bool is_in_x_range = (point[0] >= start_x && point[0] <= end_x);
    const bool is_in_y_range = (point[1] >= start_y && point[1] <= end_y);
    const bool is_in_z_range = (point[2] >= start_z && point[2] <= end_z);

    if (is_in_x_range) {
        for (Real y = start_y; y <= end_y; y += dy) {
            for (Real z = start_z; z <= end_z; z += dz) {
                point_coord.push_back(point[0]);
                point_coord.push_back(y);
                point_coord.push_back(z);
                probePoint(velocity_data, pressure_data, point[0], y, z);
            }
        }
    }

    if (is_in_y_range) {
        for (Real x = start_x; x <= end_x; x += dx) {
            for (Real z = start_z; z <= end_z; z += dz) {
                point_coord.push_back(x);
                point_coord.push_back(point[1]);
                point_coord.push_back(z);
                probePoint(velocity_data, pressure_data, x, point[1], z);
            }
        }
    }

    if (is_in_z_range) {
        for (Real x = start_x; x <= end_x; x += dx) {
            for (Real y = start_y; y <= end_y; y += dy) {
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

    //                   domain origin  +   subdomain start       + optional padding for subdomain intersection discrimination
    const Real start_x = params.originX + real(params.st_nX) * dx;
    const Real start_y = params.originY + real(params.st_nY) * dy + (!params.isOnBottom * dy);
    const Real start_z = params.originZ + real(params.st_nZ) * dz + (!params.isOnLeft * dz);

    //                  domain origin +     subdomain start     + subdomain size
    const Real end_x = params.originX + real(params.st_nX) * dx + real(nx) * dx;
    const Real end_y = params.originY + real(params.st_nY) * dy + real(ny) * dy;
    const Real end_z = params.originZ + real(params.st_nZ) * dz + real(nz) * dz;

    const bool is_in_x_range = (point[0] >= start_x && point[0] <= end_x);
    const bool is_in_y_range = (point[1] >= start_y && point[1] <= end_y);
    const bool is_in_z_range = (point[2] >= start_z && point[2] <= end_z);

    // If ask for x-axis then the point should be in range of y and z
    if (axis == 0 && is_in_y_range && is_in_z_range) {
        for (Real x = start_x; x <= end_x; x += dx) {
            point_coord.push_back(x);
            point_coord.push_back(point[1]);
            point_coord.push_back(point[2]);
            probePoint(velocity_data, pressure_data, x, point[1], point[2]);
        }
    }

    // If ask for y-axis then the point should be in range of x and z
    if (axis == 1 && is_in_x_range && is_in_z_range) {
        for (Real y = start_y; y <= end_y; y += dy) {
            point_coord.push_back(point[0]);
            point_coord.push_back(y);
            point_coord.push_back(point[2]);
            probePoint(velocity_data, pressure_data, point[0], y, point[2]);
        }
    }

    // If ask for z-axis then the point should be in range of x and y
    if (axis == 2 && is_in_x_range && is_in_y_range) {
        for (Real z = start_z; z <= end_z; z += dz) {
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
        ss << " " << pressure_values[i] << "\n";
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
