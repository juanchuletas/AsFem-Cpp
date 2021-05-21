#include<iostream>
#include "../Tools/searchelements.hpp"
#include "Asfpm.hpp"



int main(){
    //******* ONLY FOR DEVELOPERS **********
    grid_tools::froese_fischer::rmf = 5.0;
    grid_tools::froese_fischer::hmf = 1.0/24.0;
    //********* **********
    int order = 2; 
    int Ne;
    double rN = 250.0;
    std::string atom  = "He";
    int charge = 0;
    int angular = 0;
    double rc = 0.5;
    double wall = 0.0;
    std::string atomicModel = "Soft-Walls";
    std::string gridName = "Froese-Fischer";
    std::string femModel = "Fixed-Points";
    std::string integrals = "Analitic";
    std::string confType = "Free";
    int atomicN = getAtomicNumber(atom);
    double totpoints = grid_tools::froese_fischer::inverseKernel(rN,atomicN);
    int points = floor(totpoints);
    if(points%2==0){
        points = points - 1;
    }
    Ne = (points-1)/order;
    
    ASFPM asfpm{Ne,order,gridName,atomicN,charge,angular,rN};

    asfpm.setData(atomicModel,rc, wall);
    asfpm.startProgram();

    printf("****************** FEM FIXED POINTS MODEL ************************\n");
    std::cout<<"Elements: "<<Ne<<std::endl;
    std::cout<<"Points: "<<Ne*order+1<<std::endl;
    std::cout<<"Atom: "<<atom<<std::endl;
    std::cout<<"Confinement Model: "<<atomicModel<<std::endl;
    std::cout<<"Confinement Radius: "<<rc<<std::endl;
    std::cout<<"Confinement Potential: "<<wall<<std::endl;
    std::cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    printf("******************************************************************\n");

    return 0;
}