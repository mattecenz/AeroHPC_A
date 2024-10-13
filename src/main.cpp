#include "../include/utils.hpp"



int main(){
    float spatial_discr = 0.1;
    Model<float> model(1.0, 1.0, 1.0, spatial_discr);
    float nodes = model.get_nodes_number();
    //size_t nodi = model.get_nodes_number();
    std::cout << "Number of nodes: " << "dummy" << std::endl;
    std::vector<float> u0 = {1.0, 1.0, 1.0}, u1 = {2.0, 2.0, 2.0};
    std::vector<size_t> node = {1,1,1};
    model.assign_cost_field(u0);
    model.assign_cost_BC(u1);
    std::cout << "Size of u_x: " << model.u_x.size()*model.u_x[0].size() << std::endl;
    std::cout << "value ux node (2,3,4): " << model.u_x[3][21] << std::endl;
    std::cout << "du_dx: " << model.du_dx(node, spatial_discr)[0] << std::endl;
    std::cout << "du_dy: " << model.du_dy(node, spatial_discr)[0] << std::endl;
    std::cout << "du_dz: " << model.du_dz(node, spatial_discr)[0] << std::endl;
    std::cout << "laplacian u: " << model.laplacian_u(node, spatial_discr)[0] << std::endl;
    //std::cout << "p_derivative: " << model.p_derivative(node, spatial_discr) << std::endl;
    std::cout << "u_interp_x: " << model.u_interp_x(node)[0] << std::endl;
    std::cout << "u_interp_y: " << model.u_interp_y(node)[0] << std::endl;
    return 0;
}