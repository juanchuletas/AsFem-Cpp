#include<iostream>
#include<cmath>
#include "../ASFPM/Asfpm.hpp"
#include "../UserInput/UserInputData.hpp"
#include "../Tools/searchelements.hpp"
#define MIN_REQUIRED 2



int setNumberOfElements(double rN, int atomicN,int order){
    double totpoints = grid_tools::froese_fischer::inverseKernel(rN,atomicN);
    int points = floor(totpoints);
    if(points%2==0){
        points = points - 1;
    }
    int Ne = (points-1)/order;

    return Ne;
}
int help() {
   printf("Usage: asfem.exe [-s ] [-n ] [-true]\n");
   printf("\t-s: a string a\n");
   printf("\t-n: a number\n");
   printf("\t-true: a single parameter\n");

   return 1;
}

int main(int argc, char **argv){

    std::string atomSym,meshType,atomicModel,femModel,confType, integrals,multi;
    int angular,Ne,order,charge;
    double r0, rN,lambda, rC,wallValue;
    grid_tools::froese_fischer::rmf = 5.0;
    grid_tools::froese_fischer::hmf = 1.0/24.0;
    if(argc<MIN_REQUIRED){
        std::cout<<"IDIOT, YOU NEED AN INPUT FILE"<<std::endl;
        return help();
    }
    else{
        int flag = read_data(argv,femModel,Ne, order, atomicModel,confType,rC,wallValue,lambda,charge,atomSym,multi,angular,r0,rN,meshType,integrals);
        int atomicN = getAtomicNumber(atomSym);
        Ne = setNumberOfElements(rN,atomicN, order); 

        ASFPM asfpm{Ne,order,meshType,atomicN,charge,angular,rN};
        if(atomicModel=="Free" || atomicModel=="free"){
            asfpm.setData(atomicModel);
        }
        else{
            asfpm.setData(atomicModel,rC, wallValue);
            
        }
        asfpm.load(); //Assamble all the neccesary fem matrices
        asfpm.singleDiagonalization(); 
        //double *orb = asfpm.getOrbitals(allOccupied);     
        /*
            asfpm.singleDiag(); 
            double *wfn = asfpm.getWfn(full);
            double *eigenValue = asfpm.getEigenvalues(10); 

            getDensity(wfn,dens,occOrb,fenodes);
            computeHartreePotential(dens); 
            asfpm.addElectronElecron(hpot);
            asfpm.getExchangeDens(wfn); 
            if(multiplicity==singlet){
                electronicStructure::closedShell::solveHartreeFock(fij);  
            }
            else{
                electronicStructure::openShell::solveHartreeFock();
            }
        
        */
        
      
        
        printUsrData(femModel,Ne, order, atomicModel,confType,rC,wallValue,lambda,charge,atomSym,multi,angular,r0,rN,meshType,integrals);
       
    }
    return 0;
}