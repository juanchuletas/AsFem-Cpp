#include<iostream>
#include<cmath>
//#include "../ASFEM_FIXED_POINTS/asfemfp.hpp"
#include "../ASFPM/Asfpm.hpp"
#include "../UserInput/UserInputData.hpp"
#include "../Tools/searchelements.hpp"
//#include "../Static/grid_tools.hpp"
#define MIN_REQUIRED 2




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
        std::cout<<"YOU NEED AN INPUT FILE"<<std::endl;
        return help();
    }
    else{
        int flag = read_data(argv,femModel,Ne, order, atomicModel,confType,rC,wallValue,lambda,charge,atomSym,multi,angular,r0,rN,meshType,integrals);
        int atomicN = getAtomicNumber(atomSym);
        double totpoints = grid_tools::froese_fischer::inverseKernel(rN,atomicN);
        int points = floor(totpoints);
        if(points%2==0){
            points = points - 1;
        }
        Ne = (points-1)/order;
        
        // HERE COMES THE CONSTRUCTOR

        // Asfemfp myAsfem{};
        /*
        
        
        */
        
        printUsrData(femModel,Ne, order, atomicModel,confType,rC,wallValue,lambda,charge,atomSym,multi,angular,r0,rN,meshType,integrals);
       
    }
    return 0;
}