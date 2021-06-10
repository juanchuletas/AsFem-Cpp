#include<iostream>
#include<string>
#include "UserInputData.hpp"
int main (int argc, char **argv){

    std::string atomSym,meshType,atomicModel,femModel,confType, integrals;
    int angular,Ne,order,charge;
    double r0, rN,lambda, rC,wallValue;
    int flag = read_data(argv,femModel,Ne, order, atomicModel,confType,rC,wallValue,lambda,charge,atomSym,angular,r0,rN,meshType,integrals);
    std::cout<<"Your initial data:"<<std::endl;
    std::cout<<"FEM Model: "<<femModel<<std::endl;
    if(femModel=="Fixed-Elements"){
        std::cout<<"Number of Elements: "<<Ne<<std::endl;
    }
    std::cout<<"Order: "<<order<<std::endl;
    std::cout<<"Atomic Model: "<<atomicModel<<std::endl;
    if(atomicModel=="Confined"){
        if(confType=="Soft-Walls"){
            std::cout<<"Confinement: "<<confType<<std::endl;
            std::cout<<"Wall Value: "<<wallValue<<std::endl;
            std::cout<<"Rc: "<<rC<<std::endl;
        }
        else{
            std::cout<<"Confinement: "<<confType<<std::endl;
            std::cout<<"Rc: "<<rC<<std::endl;
        }
        
    }
    else if(atomicModel=="Plasma"){
        std::cout<<"Lambda Value: "<<lambda<<std::endl;
        
    }
    std::cout<<"Atom: "<<atomSym<<std::endl;
    std::cout<<"Charge: "<<charge<<std::endl;
    std::cout<<"Angular Momentum: "<<angular<<std::endl;
    std::cout<<"rInfty: "<<rN<<std::endl;

    std::cout<<"Mesh: "<<meshType<<std::endl;






}