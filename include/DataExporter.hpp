#ifndef DATAEXPORTER_HPP
#define DATAEXPORTER_HPP

#include "MPITraits.hpp"
#include "Traits.hpp"
#include <mpi.h>
#include <string>
#include <vector>

void extractData(const GridData &model,
                 std::vector<Real> &point_coord,
                 std::vector<Real> &velocity_data,
                 std::vector<Real> &pressure_data) {
    Real dx = model.structure.dx;
    Real dy = model.structure.dy;
    Real dz = model.structure.dz;

    index_t nx = model.structure.nx;
    index_t ny = model.structure.ny;
    index_t nz = model.structure.nz;

    Real origin_x = model.structure.px * dx;
    Real origin_y = model.structure.py * dy;
    Real origin_z = model.structure.pz * dz;


    for (index_t i = 0; i < nx; i++) {
        for (index_t j = 0; j < ny; j++) {
            for (index_t k = 0; k < nz; k++) {
                point_coord.push_back(origin_x + i * dx);
                point_coord.push_back(origin_y + j * dy);
                point_coord.push_back(origin_z + k * dz);
                velocity_data.push_back(model.U(i, j, k));
                velocity_data.push_back(model.V(i, j, k));
                velocity_data.push_back(model.W(i, j, k));
                pressure_data.push_back(model.P(i, j, k));
            }
        }
    }
}


template<typename T>
void swap_endian(T &var) {
    char *varArray = reinterpret_cast<char *>(&var);
    for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++)
        std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
}

void writeVtkFile(
    const std::string &filename,
    const std::string &description,
    std::vector<Real> &point_values,
    std::vector<Real> &velocity_values,
    std::vector<Real> &pressure_values)
{
    const int local_n = static_cast<int>(point_values.size() / 3);

    for (float &point_value: point_values)
        swap_endian(point_value);

    for (float &velocity_value: velocity_values)
        swap_endian(velocity_value);

    for (float &pressure_value: pressure_values)
        swap_endian(pressure_value);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_File file;
    MPI_Status status;

    const unsigned long datasize = sizeof(Real);
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

#endif //DATAEXPORTER_HPP
