#ifndef AEROHPC_A_RUNGEKUTTA_H
#define AEROHPC_A_RUNGEKUTTA_H

#include <StaggeredGrid.hpp>
#include "utils.hpp" 

using namespace utils;


#define alpha1 (64.0/120.0)
#define alpha2 (34.0/120.0)
#define alpha3 (50.0/120.0)
#define alpha4 (90.0/120.0)

//RHS function
template <typename T, Addressing_T A>
Real rhs(const StaggeredGrid<T, A> &grid, Component c, int i, int j, int k, double Kappa, int Re) {

    return Kappa*(-(conv(grid,c,i,j,k))+ (1/Re)* lap(grid,c, i,j,k));
}

//Runge-Kutta method
void rungeKutta(int n,  int Re, const StaggeredGrid<double,Addressing_T::STANDARD> grid,
StaggeredGrid<double,Addressing_T::STANDARD>  &grid_out, double deltat) {

    //Y1
    double Kappa;
    double Kappa2;
    // no default constructor.....
    StaggeredGrid<double,Addressing_T::STANDARD> K2 = grid;
    StaggeredGrid<double,Addressing_T::STANDARD> K3 = grid;

    Kappa = (alpha1*deltat) ;    

    //Y2
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                
                K2(Component::U, i,j,k) = grid(Component::U, i,j,k) + rhs(grid,Component::U,i,j,k,Kappa,Re);
                K2(Component::V, i,j,k) = grid(Component::V, i,j,k) + rhs(grid,Component::V,i,j,k,Kappa,Re);
                K2(Component::W, i,j,k) = grid(Component::W, i,j,k) + rhs(grid,Component::W,i,j,k,Kappa,Re);                
            }
        }
    }

    Kappa = (alpha3*deltat);
    Kappa2 =(alpha2*deltat);

   //Y3
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {

                K3(Component::U, i,j,k) = K2(Component::U, i,j,k) + rhs(K2,Component::U,i,j,k,Kappa,Re) - rhs(K2,Component::U,i,j,k,Kappa2,Re);
                K3(Component::V, i,j,k) = K2(Component::V, i,j,k) + rhs(K2,Component::V,i,j,k,Kappa,Re) - rhs(K2,Component::V,i,j,k,Kappa2,Re);
                K3(Component::W, i,j,k) = K2(Component::W, i,j,k) + rhs(K2,Component::W,i,j,k,Kappa,Re) - rhs(K2,Component::W,i,j,k,Kappa2,Re);  

                 
            }
        }
    }

    Kappa2 = (alpha3*deltat);
    Kappa =(alpha4*deltat);

   //Final answer
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {

                grid_out(Component::U, i,j,k) = K3(Component::U, i,j,k) + rhs(K3,Component::U,i,j,k,Kappa,Re) - rhs(grid,Component::U,i,j,k,Kappa2,Re);
                grid_out(Component::V, i,j,k) = K3(Component::V, i,j,k) + rhs(K3,Component::V,i,j,k,Kappa,Re) - rhs(grid,Component::V,i,j,k,Kappa2,Re);
                grid_out(Component::W, i,j,k) = K3(Component::W, i,j,k) + rhs(K3,Component::W,i,j,k,Kappa,Re) - rhs(grid,Component::W,i,j,k,Kappa2,Re);  

                 
            }
        }
    }


}

   

#endif //AEROHPC_A_RUNGEKUTTA_H
