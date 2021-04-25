#include<iostream>
#include "grid_tools.hpp"






int main (){

    grid_tools::froese_fischer::hmf = 1.0/32.0;
    grid_tools::froese_fischer::rmf = 5.0;
    int order = 2;
    int nele;

    double value = grid_tools::froese_fischer::kernel(383,2.0);

    double totpoints = grid_tools::froese_fischer::inverseKernel(20.0,2.0);


    int points = floor(totpoints);
    if(points%2==0){
        points = points - 1;
    }
    nele = (points-1)/order;

    std::cout<<"Ne is : "<<nele<<std::endl;



    return 0;
}