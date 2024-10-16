#include <utils.hpp>
#include <iostream>

using namespace std;

int main(int argc, char** argv){
    //Random stuff to see if everything works (hopefully correctly)
    constexpr int lx=3;
    constexpr int ly=3;
    constexpr int lz=3;

    StaggeredGrid<Real,Addressing_T::STANDARD> grid(lx,ly,lz);
    for(int i=0;i<lx;i++){
        for(int j=0;j<ly;j++){
            for(int k=0;k<lz;k++){
                grid(Component::U,i,j,k) = i*j*k;
                grid(Component::V,i,j,k) = i*j*k;
                grid(Component::W,i,j,k) = i*j*k;
                grid(Component::P,i,j,k) = i*j*k;
            }
        }
    }


    cout << utils::d_dx<Real,Addressing_T::STANDARD>(grid,Component::U,1,1,1) << endl;
    cout << utils::d2_dz2<Real,Addressing_T::STANDARD>(grid,Component::W,1,1,1) << endl;
    cout << utils::get_interpolation<Real,Addressing_T::STANDARD>(grid,Component::U,Component::V,1,1,1) << endl;
    return 0;
}