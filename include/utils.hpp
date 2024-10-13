#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>

template<typename T>
class Model{
    public:
    Model(const T x_length, const T y_length, const T z_length, const T spatial_discr):
        x_Length(x_length), y_Length(y_length), z_Length(z_length), spatial_discr(spatial_discr)
    {
       size_t nodes = get_nodes_number();

    };

    size_t get_nodes_number();

    std::vector<T> du_dx(std::vector<size_t>& node, T spatial_discr);
    std::vector<T> du_dy(std::vector<size_t>& node, T spatial_discr);
    std::vector<T> du_dz(std::vector<size_t>& node, T spatial_discr);
    std::vector<T> p_derivative(std::vector<size_t>& node, T spatial_discr);

    std::vector<T> u_interp_x(std::vector<size_t>& node);
    std::vector<T> u_interp_y(std::vector<size_t>& node);
    std::vector<T> u_interp_z(std::vector<size_t>& node);

    std::vector<T> laplacian_u(std::vector<size_t>& node, T spatial_discr);

    void assign_cost_field(std::vector<T> u0);

    void assign_cost_BC(std::vector<T> u0);



    std::vector<std::vector<T>> u_x, u_y, u_z, pressure;
    size_t x_nodes, y_nodes, z_nodes;

    private:
    T x_Length, y_Length, z_Length, spatial_discr;
    

};






#endif //UTILS_HPP