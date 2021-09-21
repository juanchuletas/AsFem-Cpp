#include<iostream>
#include "Asfem.hpp"

//***************** ASFEM CONSTRUCTORS & DESTRUCTORS ****************************************************************
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
    vr_pot = nullptr;
    vh_pot = nullptr;
}
template<class T>
Asfem<T>::Asfem(double inr0, double inrN, int inNe,int inOrder,std::string inatomName,std::string nameMesh,std::string inputModel,std::string intAtomicModel,std::string potInt)
:r0{inr0},rN{inrN},Ne{inNe},order{inOrder},meshType{nameMesh},femModel{inputModel},atomName{inatomName},atomicModel{intAtomicModel},integrationScheme{potInt}{
    std::cout<<"AFEM CONSTRUCTOR\n";
    polynomial = order + 1;
    totalNodes = Ne*order + 1;
    boundaryNodes = totalNodes-2;
    atomicN = getAtomicNumber(atomName);
    femGrid.setGridData(r0,rN,Ne,order,meshType,atomicN);
    femGrid.createGrid();
    femStuff.setFemData(Ne,order);
    globMatSize = totalNodes*totalNodes;
    sij.setMatrix(boundaryNodes*boundaryNodes);
    kij.setMatrix(boundaryNodes*boundaryNodes);
    lij.setMatrix(boundaryNodes*boundaryNodes);
    vr_pot = new T[boundaryNodes];
    vh_pot = new T[boundaryNodes];

}
template<class T>
Asfem<T>::Asfem(double inr0, double inrN,int inOrder,std::string inatomName,std::string nameMesh,std::string inputModel,std::string intAtomicModel,std::string potInt)
:r0{inr0},rN{inrN},order{inOrder},meshType{nameMesh},femModel{inputModel},atomName{inatomName},atomicModel{intAtomicModel},integrationScheme{potInt}{
    atomicN = getAtomicNumber(atomName);
    Ne = computeFixedElements();
    polynomial = order + 1;
    totalNodes = Ne*order + 1;
    boundaryNodes = totalNodes-2;
    femGrid.setGridData(r0,rN,Ne,order,meshType,atomicN);
    femGrid.createGrid();
    femStuff.setFemData(Ne,order);
    globMatSize = totalNodes*totalNodes;
    sij.setMatrix(boundaryNodes*boundaryNodes);
    kij.setMatrix(boundaryNodes*boundaryNodes);
    lij.setMatrix(boundaryNodes*boundaryNodes);
  
}
template<class T>
Asfem<T>::~Asfem(){

    std::cout<<"Destructor\n";
    delete [] vh_pot;
    delete [] vr_pot;
}
//***************** END OF ASFEM CONSTRUCTORS & DESTRUCTORS ****************************************************************
//
//****************** ASFEM PRIVATE METHODS *********************************************************************************
template<class T>
void Asfem<T>::applyBoundaryConditions(T* s_matG,T* k_matG,T* l_matG){
    int l=0;
    int n = totalNodes;
    for(int i=1;i<n-1;i++)
    {
        for(int j=1;j<n-1;j++)
        {
            sij[l] = s_matG[i*n+j]; //Overlap matrix
            kij[l] = k_matG[i*n+j]; //Kinect Matrix
            lij[l] = l_matG[i*n+j]; //poisson matrix
            l++;
        }
    }
    sij.printMatrix(boundaryNodes,boundaryNodes,"Overlap Matrix");
    kij.printMatrix(boundaryNodes,boundaryNodes,"Kinect Matrix");
    lij.printMatrix(boundaryNodes,boundaryNodes,"Poisson's Matrix");

}
template<class T>
void Asfem<T>::integrateFembasis(){
//****** GLOBAL MATRICES *************
T *s_matG = new T[totalNodes*totalNodes];
T *k_matG = new T[totalNodes*totalNodes];
T *l_matG = new T[totalNodes*totalNodes];
//**** ELEMENTAL MATRICES ***************
T *feMatS = new T[polynomial*polynomial];
T *feMatK = new T[polynomial*polynomial];
T *feMatL = new T[polynomial*polynomial];
//**************************************

feMatS = femStuff.getOverlap();
feMatK = femStuff.getKinect();
feMatL = femStuff.getKinect();
double coeff,coeff2;

if(femModel=="Fixed Elements"){
    feMatS = femStuff.getOverlap();
    feMatK = femStuff.getKinect();
    feMatL = femStuff.getKinect();
    double coeff,coeff2;  
    for(int ei=0; ei<Ne; ei++)
    {
        coeff = 0.5*femGrid.getElementSize(ei);
        coeff2 = femGrid.getElementSize(ei);
        for(int nu=0; nu<polynomial; nu++)
        {
            int index_nu = polynomial*ei + nu;
            int l = femStuff.getLinkMatIndex(index_nu);
            for(int mu = 0; mu<polynomial; mu++)
            {
                int index_mu = polynomial*ei+mu;
                int m = femStuff.getLinkMatIndex(index_mu);
                s_matG[l*totalNodes+m] += coeff*feMatS[polynomial*nu+mu];
                k_matG[l*totalNodes+m] += feMatK[polynomial*nu+mu]/coeff;
                l_matG[l*totalNodes+m] += 2.0*feMatL[polynomial*nu+mu]/coeff2;
            }
        }
    }
}
else
{
    double x[polynomial];
    double v[polynomial];
    for(int ei=0; ei<Ne; ei++)
    {
        for(int j=0; j<polynomial; j++)
        {
            int indx = polynomial*ei + j;
            int i = femStuff.getLinkMatIndex(indx);
            x[j] = femGrid.getGridItem(i);
        }
        feMatS = femStuff.overlapMod(x);
        feMatK = femStuff.kinectMod(x);
        for(int nu=0; nu<polynomial; nu++)
        {
            int index_nu = polynomial*ei + nu;
            int l = femStuff.getLinkMatIndex(index_nu);
            for(int mu=0; mu<polynomial; mu++)
            {
                int index_nu = polynomial*ei+mu;
                int m = femStuff.getLinkMatIndex(index_mu);
                s_matG[l*totalNodes+m] += feMatS[polynomial*nu+mu];
                k_matG[l*totalNodes+m] += feMatK[polynomial*nu+mu];
            }
        }
    }
}

applyBoundaryConditions(s_matG,k_matG,l_matG);

delete [] feMatS;
delete [] feMatK;
delete [] feMatL;
delete [] s_matG;
delete [] k_matG;
delete [] l_matG;
}
template<class T>
int Asfem<T>::computeFixedElements(){
    int nele;
    if(meshType=="Froese-Fischer"){
        double totpoints = atomicN*32.0*(log(rN) + 5);  
        int points = floor(totpoints);
        if(points%2==0){
            points = points - 1;
        }
        nele = (points-1)/order;
    }

    return nele;
    
}
template<class T>
void Asfem<T>::evaluateExternalPotential(){
    if(atomicModel=="Free-Atom"){
        for(int i=0; i<Ne; i++){
            vr_pot[i] = externalPotential(femGrid[i]);
            printf("V(r) = %lf\n",vr_pot[i]);
        }
    }
}
//****************** END OF ASFEM PRIVATE METHODS *********************************************************************************
// *****************
// ****************  ASFEM PUBLIC METHODS ******************************************************************************************
template<class T>
void Asfem<T>::runProgram(){

    //femGrid.createGrid();
    //femGrid.printGrid();
    integrateFembasis();
    

}

template<class T>
void Asfem<T>::printInputData(){

    std::cout<<"********** ASFEM INITIAL DATA **********\n";
    std::cout<<"Finite Element Model: "<<femModel<<std::endl;
    std::cout<<"Number of Elements: "<<Ne<<std::endl;
    std::cout<<"Total Points: "<<(Ne*order + 1)<<std::endl;
    std::cout<<"Polynomial Order: "<<order<<std::endl;
    std::cout<<"Atom: "<<atomName<<std::endl;
    std::cout<<"Model: "<<atomicModel<<std::endl;
    if(atomicModel=="Confined-Atom"){
        std::cout<<"Confinement: "<<confType<<std::endl;
        if(confType=="Soft-Walls"){
            std::cout<<"Cut Radii: "<<cutRad<<std::endl;
            std::cout<<"Wall-Value: "<<wallValue<<std::endl;
        }
        else{
            std::cout<<"Cut Radii: "<<cutRad<<std::endl;
        }
    }
    else if(atomicModel=="Plasma"){
        std::cout<<"Lambda Value: "<<lambda<<std::endl;
    }
    std::cout<<"Practical Infinity: "<<rN<<std::endl;
    std::cout<<"Grid Form: "<<meshType<<std::endl;
    std::cout<<"External Potential Integral: "<<integrationScheme<<std::endl;
    std::cout<<"***************************************\n";
}
// ****************  END OF ASFEM PUBLIC METHODS ******************************************************************************************
