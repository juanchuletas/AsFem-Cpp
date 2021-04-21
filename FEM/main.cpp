#include <iostream>
#include<string>
#include "FiniteElement.hpp"




int main(){
    int Ne = 10; 
    int order = 2;
    double *sij,*kij;
    sij = new double[(order+1)*(order+1)];
    kij = new double[(order+1)*(order+1)];
    std::string name = "Froese-Fischer";
    FEM fem;
    fem.setFemData(Ne,order);
    /* sij = fem.getOverlap();
    kij = fem.getKinect();
    printf("Overlap[0] = %lf\n",sij[2]);
    for(int i=0; i<(order+1)*(order+1); i++){ 
        printf("%lf\n", kij[i]);
    }
    /* for(int i=0; i<Ne*order+1; i++){

        std::cout<<fem.getLinkMatIndex(i)<<std::endl;
    } */ 



    delete [] sij;
    delete [] kij;

}