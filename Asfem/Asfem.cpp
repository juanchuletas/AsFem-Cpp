#include<iostream>
#include "Asfem.hpp"


template<class T>
Asfem<T>::Asfem(){
    Ne = 0;
    order = 0;
    totalNodes = 0;
    polynomial = 0;
    boundaryNodes = 0;
    atomicN = 0;
    globMatSize = 0;
    r0=0.0;
    rN = 0.0;
}
template<class T>
Asfem<T>::Asfem(double inr0, double inrN, int inNe,int inOrder,int atom,std::string nameMesh)
:r0{inr0},rN{inrN},Ne{inNe},order{inOrder},meshType{nameMesh},atomicN{atom}{
    std::cout<<"AsFem CONSTRUCTOR\n";
    polynomial = order + 1;
    totalNodes = Ne*order + 1;
    boundaryNodes = totalNodes-2;
    femGrid.setGridData(r0,rN,Ne,order,meshType,atomicN);
    femGrid.createGrid();
    femStuff.setFemData(Ne,order);
    globMatSize = totalNodes*totalNodes;
    s_mat.setMatrix(globMatSize);
  
}
template<class T>
Asfem<T>::~Asfem(){

    std::cout<<"Destructor\n";

}
template<class T>
void Asfem<T>::integrateFembasis(){
T *feMatS = new T[polynomial*polynomial];
feMatS = femStuff.getOverlap();
double coeff;
for(int ei=0; ei<Ne; ei++)
{
    coeff = 0.5*femGrid.getElementSize(ei);
    printf("%lf\n",femGrid[1]);
    
    for(int nu=0; nu<polynomial; nu++){
        int index_nu = polynomial*ei + nu;
        int l = femStuff.getLinkMatIndex(index_nu);
        for(int mu = 0; mu<polynomial; mu++){

            int index_nu = polynomial*ei+nu;
            int m = femStuff.getLinkMatIndex(index_nu);
            s_mat[l*totalNodes+m] += coeff*feMatS[polynomial*nu+mu];


        }
    }
}

s_mat.printMatrix(totalNodes,totalNodes,"overlap matrix");

/* T *feMatK = T[polynomial*polynomial];
T *feMatP = T[polynomial*polynomial];
 */

delete [] feMatS;

}
template<class T>
void Asfem<T>::assambleMatrices(){

    integrateFembasis();
}
template<class T>
void Asfem<T>::runProgram(){

    femGrid.createGrid();
    assambleMatrices();

}