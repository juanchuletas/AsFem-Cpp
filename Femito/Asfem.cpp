#include<iostream>
#include "Asfem.hpp"

template<class T>
Asfem<T>::Asfem(double inr0, double inrN, int inNe,int inOrder,int atom,std::string nameMesh)
:r0{inr0},rN{inrN},Ne{inNe},order{inOrder},meshType{nameMesh},atomicN{atom}{
    std::cout<<"FEMITO CONSTRUCTOR\n";
    polynomial = order + 1;
    totalNodes = Ne*order + 1;
    boundaryNodes = totalNodes-2;
    //femGrid.createGrid(meshType,atomicN);
    femStuff.setFemData(Ne,order,meshType,atomicN);
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
double coeff;
for(int ei=0; ei<Ne; ei++)
{
    coeff = 0.5*femStuff.getElementSize(ei);
    
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


}
template<class T>
void Asfem<T>::runProgram(){

    assambleMatrices();


}