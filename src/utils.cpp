#include "../include/utils.hpp" 


// Defining mesh parameters****************************************************************************************

/**
 * vectors are two components: the first indicate the z direction, the second is a row-major vector that stores the x/y components
*/

template<typename T>
size_t Model<T>::get_nodes_number() {
         x_nodes = static_cast<size_t>(floor(x_Length/spatial_discr));
         y_nodes = static_cast<size_t>(floor(y_Length/spatial_discr));
         z_nodes = static_cast<size_t>(floor(z_Length/spatial_discr));

        size_t nodes = x_nodes*y_nodes*z_nodes;

        u_x.resize(z_nodes);
        u_y.resize(z_nodes);
        u_z.resize(z_nodes);
        pressure.resize(z_nodes);

        for(size_t i = 0; i < z_nodes; i++){
            u_x[i].resize(y_nodes*x_nodes);
            u_y[i].resize(y_nodes*x_nodes);
            u_z[i].resize(y_nodes*x_nodes);
            pressure[i].resize(y_nodes*x_nodes);
        }

        
        

        return x_nodes*y_nodes*z_nodes;
    }

    template size_t Model<float>::get_nodes_number();
    template size_t Model<double>::get_nodes_number();

    //Evaluate the derivative of the velocity****************************************************************************************

    /**
     * with the "node" parameter i mean a 3 component vector that define a single node in position (x,y,z), onche defined the node,  
     * the function evaluate the derivative of the velocity in the x direction, y direction and z direction, the result is stored 
     * in a vector of 3 components 
    */

    template<typename T>
    std::vector<T> Model<T>::du_dx(std::vector<size_t>& node, T spatial_discr){
        std::vector<T> result;
        result.resize(3);
        //dux_dx
        result[0] = (u_x[node[2]][node[1]*x_nodes+node[0]+1] - u_x[node[2]][node[1]*x_nodes+node[0]-1])/(2*spatial_discr);
        //duy_dx
        result[1] = (u_y[node[2]][(node[1])*x_nodes+node[0]+1] - u_y[node[2]][(node[1])*x_nodes+node[0]-1])/(2*spatial_discr);
        //duz_dx
        result[2] = (u_z[node[2]][(node[1])*x_nodes+node[0]+1] - u_z[node[2]][(node[1])*x_nodes+node[0]-1])/(2*spatial_discr);

        return result;
        
    }

    template std::vector<float> Model<float>::du_dx(std::vector<size_t>& node, float spatial_discr);
    template std::vector<double> Model<double>::du_dx(std::vector<size_t>& node, double spatial_discr);

    template<typename T>
    std::vector<T> Model<T>::du_dy(std::vector<size_t>& node, T spatial_discr){
        std::vector<T> result;
        result.resize(3);
        //dux_dy
        result[0] = (u_x[node[2]][(node[1]+1)*x_nodes+node[0]] - u_x[node[2]][(node[1]-1)*x_nodes+node[0]])/(2*spatial_discr);
        //duy_dy
        result[1] = (u_y[node[2]][(node[1]+1)*x_nodes+node[0]] - u_y[node[2]][(node[1]-1)*x_nodes+node[0]])/(2*spatial_discr);
        //duz_dy
        result[2] = (u_z[node[2]][(node[1]+1)*x_nodes+node[0]] - u_z[node[2]][(node[1]-1)*x_nodes+node[0]])/(2*spatial_discr);

        return result;
        
    }

    template std::vector<float> Model<float>::du_dy(std::vector<size_t>& node, float spatial_discr);
    template std::vector<double> Model<double>::du_dy(std::vector<size_t>& node, double spatial_discr);

    template<typename T>
    std::vector<T> Model<T>::du_dz(std::vector<size_t>& node, T spatial_discr){
        std::vector<T> result;
        result.resize(3);
        //dux_dz
        result[0] = (u_x[node[2]+1][(node[1])*x_nodes+node[0]] - u_x[node[2]-1][(node[1])*x_nodes+node[0]])/(2*spatial_discr);
        //duy_dz
        result[1] = (u_y[node[2]+1][(node[1])*x_nodes+node[0]] - u_y[node[2]-1][(node[1])*x_nodes+node[0]])/(2*spatial_discr);
        //duz_dz
        result[2] = (u_z[node[2]+1][(node[1])*x_nodes+node[0]] - u_z[node[2]-1][(node[1])*x_nodes+node[0]])/(2*spatial_discr);

        return result;
        
    }

    template std::vector<float> Model<float>::du_dz(std::vector<size_t>& node, float spatial_discr);
    template std::vector<double> Model<double>::du_dz(std::vector<size_t>& node, double spatial_discr);

    template<typename T>
    std::vector<T> Model<T>::p_derivative(std::vector<size_t>& node, T spatial_discr){
        std::vector<T> result;
        result.resize(3);
        //dp_dx
        result[0] = (pressure[node[2]][(node[1])*x_nodes+node[0]+1] - pressure[node[2]][(node[1])*x_nodes+node[0]-1])/(2*spatial_discr);
        //dp_dy
        result[1] = (pressure[node[2]][(node[1]+1)*x_nodes+node[0]] - pressure[node[2]][(node[1]-1)*x_nodes+node[0]])/(2*spatial_discr);
        //dp_dz
        result[2] = (pressure[node[2]+1][(node[1])*x_nodes+node[0]] - pressure[node[2]-1][(node[1])*x_nodes+node[0]])/(2*spatial_discr);

        return result;
        
    }

    template std::vector<float> Model<float>::p_derivative(std::vector<size_t>& node, float spatial_discr);
    template std::vector<double> Model<double>::p_derivative(std::vector<size_t>& node, double spatial_discr);

    //Interpolation of the velocity field****************************************************************************************

    template<typename T>
    std::vector<T> Model<T>::u_interp_x(std::vector<size_t>& node){
        std::vector<T> result;
        result.resize(2);
        //v interpolated in x
        result[0] = (u_y[node[2]][(node[1])*x_nodes+node[0]] + u_y[node[2]][(node[1]-1)*x_nodes+node[0]] + u_y[node[2]][(node[1])*x_nodes+node[0]+1] + u_y[node[2]][(node[1]-1)*x_nodes+node[0]+1])/4;
        //w interpolated in x
        result[1] = (u_z[node[2]][(node[1])*x_nodes+node[0]] + u_z[node[2]-1][(node[1])*x_nodes+node[0]] + u_z[node[2]][(node[1])*x_nodes+node[0]+1] + u_z[node[2]-1][(node[1])*x_nodes+node[0]+1])/4;

        return result;

    }

    template std::vector<float> Model<float>::u_interp_x(std::vector<size_t>& node);
    template std::vector<double> Model<double>::u_interp_x(std::vector<size_t>& node);

    template<typename T>
    std::vector<T> Model<T>::u_interp_y(std::vector<size_t>& node){
        std::vector<T> result;
        result.resize(2);
        //u interpolated in y
        result[0] = (u_x[node[2]][(node[1])*x_nodes+node[0]] + u_x[node[2]][(node[1])*x_nodes+node[0]-1] + u_x[node[2]][(node[1]+1)*x_nodes+node[0]] + u_x[node[2]][(node[1]+1)*x_nodes+node[0]-1])/4;
        //w interpolated in y
        result[1] = (u_z[node[2]][(node[1])*x_nodes+node[0]] + u_z[node[2]-1][(node[1])*x_nodes+node[0]] + u_z[node[2]][(node[1]+1)*x_nodes+node[0]] + u_z[node[2]-1][(node[1]+1)*x_nodes+node[0]])/4;

        return result;

    }

    template std::vector<float> Model<float>::u_interp_y(std::vector<size_t>& node);
    template std::vector<double> Model<double>::u_interp_y(std::vector<size_t>& node);

    template<typename T>
    std::vector<T> Model<T>::u_interp_z(std::vector<size_t>& node){
        std::vector<T> result;
        result.resize(2);
        //u interpolated in z
        result[0] = (u_x[node[2]][(node[1])*x_nodes+node[0]] + u_x[node[2]][(node[1])*x_nodes+node[0]-1] + u_x[node[2]+1][(node[1])*x_nodes+node[0]] + u_x[node[2]+1][(node[1])*x_nodes+node[0]-1])/4;
        //v interpolated in z
        result[1] = (u_y[node[2]][(node[1])*x_nodes+node[0]] + u_y[node[2]][(node[1]-1)*x_nodes+node[0]] + u_y[node[2]+1][(node[1])*x_nodes+node[0]] + u_y[node[2]+1][(node[1]-1)*x_nodes+node[0]])/4;

        return result;

    }

    template std::vector<float> Model<float>::u_interp_z(std::vector<size_t>& node);
    template std::vector<double> Model<double>::u_interp_z(std::vector<size_t>& node);

    //Evaluate the laplacian of the velocity field****************************************************************************************

    template<typename T>
    std::vector<T> Model<T>::laplacian_u(std::vector<size_t>& node, T spatial_discr){
        std::vector<T> result;
        result.resize(3);
        //laplacian u
        result[0] = (u_x[node[2]][(node[1])*x_nodes+node[0]+1] - 2*u_x[node[2]][(node[1])*x_nodes+node[0]] + u_x[node[2]][(node[1])*x_nodes+node[0]-1])/(spatial_discr*spatial_discr) + (u_x[node[2]][(node[1]+1)*x_nodes+node[0]] - 2*u_x[node[2]][(node[1])*x_nodes+node[0]] + u_x[node[2]][(node[1]-1)*x_nodes+node[0]])/(spatial_discr*spatial_discr) + (u_x[node[2]+1][(node[1])*x_nodes+node[0]] - 2*u_x[node[2]][(node[1])*x_nodes+node[0]] + u_x[node[2]-1][(node[1])*x_nodes+node[0]])/(spatial_discr*spatial_discr);
        //laplacian v
        result[1] = (u_y[node[2]][(node[1])*x_nodes+node[0]+1] - 2*u_y[node[2]][(node[1])*x_nodes+node[0]] + u_y[node[2]][(node[1])*x_nodes+node[0]-1])/(spatial_discr*spatial_discr) + (u_y[node[2]][(node[1]+1)*x_nodes+node[0]] - 2*u_y[node[2]][(node[1])*x_nodes+node[0]] + u_y[node[2]][(node[1]-1)*x_nodes+node[0]])/(spatial_discr*spatial_discr) + (u_y[node[2]+1][(node[1])*x_nodes+node[0]] - 2*u_y[node[2]][(node[1])*x_nodes+node[0]] + u_y[node[2]-1][(node[1])*x_nodes+node[0]])/(spatial_discr*spatial_discr);
        //laplacian w
        result[2] = (u_z[node[2]][(node[1])*x_nodes+node[0]+1] - 2*u_z[node[2]][(node[1])*x_nodes+node[0]] + u_z[node[2]][(node[1])*x_nodes+node[0]-1])/(spatial_discr*spatial_discr) + (u_z[node[2]][(node[1]+1)*x_nodes+node[0]] - 2*u_z[node[2]][(node[1])*x_nodes+node[0]] + u_z[node[2]][(node[1]-1)*x_nodes+node[0]])/(spatial_discr*spatial_discr) + (u_z[node[2]+1][(node[1])*x_nodes+node[0]] - 2*u_z[node[2]][(node[1])*x_nodes+node[0]] + u_z[node[2]-1][(node[1])*x_nodes+node[0]])/(spatial_discr*spatial_discr);

        return result;

    }

    template std::vector<float> Model<float>::laplacian_u(std::vector<size_t>& node, float spatial_discr);
    template std::vector<double> Model<double>::laplacian_u(std::vector<size_t>& node, double spatial_discr); 

    //Assign the boundary conditions****************************************************************************************

    template<typename T>
    void Model<T>::assign_cost_BC(std::vector<T> u0){
        for(size_t i = 0; i < u_x[0].size(); i++){
            u_x[0][i] = u0[0];
            u_y[0][i] = u0[1];
            u_z[0][i] = u0[2];
        }
    }
            

    template void Model<float>::assign_cost_BC(std::vector<float> u0);
    template void Model<double>::assign_cost_BC(std::vector<double> u0);

    //Assign the cost field****************************************************************************************

    template<typename T>
    void Model<T>::assign_cost_field(std::vector<T> u0){
        for(size_t i = 0; i < u_x[0].size(); i++){
            for(size_t j = 0; j < z_nodes; j++){
                u_x[j][i] = u0[0];
                u_y[j][i] = u0[1];
                u_z[j][i] = u0[2];
            }   
        }
    }

    template void Model<float>::assign_cost_field(std::vector<float> u0);
    template void Model<double>::assign_cost_field(std::vector<double> u0);