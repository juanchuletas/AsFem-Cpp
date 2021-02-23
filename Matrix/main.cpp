#include<iostream>
#include <iomanip>
#include "Matrix.cpp"



int main (){
    int Ne = 5;
    double value = 7.5;
    Matrix<double> grid;
    Matrix<double> elements{Ne,value};
    Matrix<double> nodes{Ne,value}; //Ten nodes array;
    //nodes.fillWithZeros();
    //nodes[0] = 3.9;
    grid.setMatrix(Ne);
    for(int i=0; i<grid.getSize(); i++){
        grid[i] = elements[i];
        std::cout<<grid[i]<<std::endl;

    }
    //printf("%lf\n",nodes[21]);


    return 0;
}