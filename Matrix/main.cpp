#include<iostream>
#include <iomanip>
#include "Matrix.hpp"



void foo(double *a,int N){

    std::cout<<"Foo\n";
    for(int i=0; i<N; i++){
        std::cout<<a[i]<<std::endl;
    }
}


int main (){
    int Ne = 5;
    double value = 7.5;
    double tryit[] = {1.4,5.6,7.5,1.9,1.1};
    Matrix<double> grid;
    Matrix<double> elements{Ne,value};
    Matrix<double> nodes{Ne,value}; //Ten nodes array;
    //nodes.fillWithZeros();
    //nodes[0] = 3.9;
    grid.setMatrix(Ne);
    for(int i=0; i<grid.getSize(); i++){
        grid[i] = tryit[i];
        std::cout<<grid[i]<<std::endl;
    }
    grid.sort();
     for(int i=0; i<grid.getSize(); i++){
        std::cout<<grid[i]<<std::endl;
    }
    foo(&grid[0],grid.getSize());
    //Matrix<double> = arr{1,2,3};
    //printf("%lf\n",nodes[-1]);



    return 0;
}