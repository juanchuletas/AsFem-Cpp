#if !defined(_FEMITO_H_)
#define _FEMITO_H_
#include<iostream>
#include <string>
#include<cmath>
#include "Elements/searchelements.hpp"
#include "../FEM/FiniteElement.hpp"
template<class T> class Asfem {
    int Ne;
    int order;
    int totalNodes;
    int polynomial;
    int boundaryNodes;
    int atomicN;
    int globMatSize;
    T r0,rN,wallValue,cutRad,lambda;
    std::string meshType;
    std::string femModel;
    std::string atomName;
    std::string atomicModel;
    std::string confType;
    std::string integrationScheme;


    Matrix<T> fij,hij,sij,kij,lij,wfn,eigen_val;
    T *vr_pot, *vh_pot;
    Grid<T> femGrid; //grid for FEM
    FEM<T> femStuff;
    private:
        //Internal Methods: Not user access    
        void integrateFembasis();
        void applyBoundaryConditions(T *,T *, T*);
        int computeFixedElements();
        void evaluateExternalPotential();
        /* void computeHartreePotential(/*Calls FEM poisson solver?);
        void electronicStructure(/*bunch of matrices to fill? );
        void HatreeFock(/*RHF or UHF ); //Hartree Fock Method
        void DFT(/*have no idea); //Density functional methods */

    public:
        Asfem();
        Asfem(double inr0, double inrN, int inNe,int inOrder,std::string atomName,std::string nameMesh,std::string inputModel,std::string intAtomicModel,std::string potInteg); //Default constructor
        Asfem(double inr0, double inrN, int inOrder,std::string atomName,std::string nameMesh,std::string inputModel,std::string intAtomicModel,std::string potInteg); //Default constructor
        ~Asfem();
        // Public Methods
        void runProgram();
        //void getUserData();
        void printFinalResults();
        void getUsrData();
        void printInputData();

        /* void HartreeFock();
        void SCF();
        void PoissonSolver(); */






};

#endif // _FEMITO_H_
