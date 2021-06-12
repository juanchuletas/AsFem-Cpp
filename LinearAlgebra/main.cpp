#include<iostream>
#include "RealSymMatrix.hpp"


int main(){
    int N = 3;
    int tot = N*N;
    double a[] = {1,1,1,1,4,3,1,3,1};

    Matrix<double> sij;
    sij.setMatrix(tot);

    for(int i = 0; i<tot; i++){
        sij[i] = a[i];
        //std::cout<<a[i]<<std::endl;
    }
    sij.printMatrix(N,N,"Dummy");

    RealSymMat sym{sij,N,N};
    sym.printResult();


    return 0;
}