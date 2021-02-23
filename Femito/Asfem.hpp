#if !defined(_FEMITO_H_)
#define _FEMITO_H_
#include<iostream>
#include <string>
#include "../FEM/FiniteElement.cpp"
template<class T> class Asfem {
    int Ne;
    int order;
    int totalNodes;
    int polynomial;
    int boundaryNodes;
    int atomicN;
    int globMatSize;
    double r0,rN;
    std::string meshType;

    Matrix<T> wfn,eigenVal,rho;
    Matrix<T> s_mat;
    Matrix<T> extPot;
    Matrix<T> vhPot;
    //Grid<T> femGrid; //grid for FEM
    FEM<T> femStuff;
    private:    
        void assambleMatrices();
        void integrateFembasis();


    public:
        Asfem(double r0, double rN, int Ne, int order,int atom,std::string meshName); //Default constructor
        ~Asfem();
        void runProgram();
};

#endif // _FEMITO_H_
