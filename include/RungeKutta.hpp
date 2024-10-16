#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include <StaggeredGrid.hpp>
#include "utils.hpp" 

using namespace utils;


#define alpha1 (64.0/120.0)
#define alpha2 (30.0/120.0)
#define alpha3 (50.0/120.0)


void compute (){

} 


//Runge-Kutta method
void rungeKutta(int n,  int Re, StaggeredGrid<double,Addressing_T::STANDARD> grid,
StaggeredGrid<double,Addressing_T::STANDARD>  &grid_out, double deltat, double h) {

    StaggeredGrid<double,Addressing_T::STANDARD> grid_out = grid;

      for (int i = 1; i < n; ++i) {
        for (int j = 1; j < n; ++j) {
            for (int k = 1; k < n; ++k) {

                double Kappa = (alpha*deltat)/Re ;

                 grid_out(Component::U, i,j,k) =  grid_out(Component::U, (i-1), (j-1), (k-1));

                 grid_out(Component::V, i,j,k) = grid_out(Component::V, (i-1), (j-1), (k-1));

                 grid_out(Component::W, i,j,k) = grid_out(Component::W, (i-1), (j-1), (k-1));
            }
        }
    }
    

}

//RHS function
template <typename T, Addressing_T A>
void rhs(const StaggeredGrid<T, A> grid, const StaggeredGrid<T, A> &grid_out,  Component c, int i, int j, int k, double Kappa) {

    
    switch (c){
        case Component::U:

            grid_out(Component::U, i,j,k) =  grid(c,i,j,k)*d_dx(grid,c,i,j,k) + 
            get_interpolation(grid, Component::U, Component::V, i, j, k)*d_dy(grid,c,i,j,k) 
            + get_interpolation(grid, Component::U, Component::W, i, j, k)*d_dz(grid,c,i,j,k) 
            + Kappa*(d2_dx2(grid, Component::U, i, j, k)+d2_dy2(grid, Component::U, i, j, k)+d2_dz2(grid, Component::U, i, j, k));

            break;

        case Component::V:
            grid_out(Component::V, i,j,k) =  get_interpolation(grid, Component::V, Component::U, i, j, k)*d_dx(grid,c,i,j,k) 
            + grid(c,i,j,k)*d_dy(grid,c,i,j,k) + 
            get_interpolation(grid, Component::V, Component::W, i, j, k)*d_dz(grid,c,i,j,k)
            +  Kappa*(d2_dx2(grid, Component::V, i, j, k)+d2_dy2(grid, Component::V, i, j, k)+d2_dz2(grid, Component::V, i, j, k));

            break;

        case Component::W:
            grid_out(Component::W, i,j,k) = get_interpolation(grid, Component::W, Component::U, i, j, k)*d_dx(grid,c,i,j,k) 
            + get_interpolation(grid, Component::W, Component::V, i, j, k)*d_dy(grid,c,i,j,k)
            + grid(c,i,j,k)*d_dz(grid,c,i,j,k)
            +  Kappa*(d2_dx2(grid, Component::W, i, j, k)+d2_dy2(grid, Component::W, i, j, k)+d2_dz2(grid, Component::W, i, j, k));

            break;

        }
}

#endif //AEROHPC_A_RUNGEKUTTA_H