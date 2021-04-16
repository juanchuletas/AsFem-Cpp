#include<iostream>
#include<string>
#include "Asfem.cpp"



int main ()
{


    int Ne;
    int order = 2;
    double r0 = 0.0;
    double rN = 20.0;
    int poissNe = 100;
    int angular = 0;
    std::string atom = "H";
    std::string atomicModel = "Free-Atom";
    std::string gridName = "Froese-Fischer";
    std::string femModel = "Fixed Points";
    std::string integrals = "Analitic";
    int atomicN = getAtomicNumber(atom);
    if(femModel=="Fixed Points"){
        printf("Fixed Points\n");
        int atomicN = getAtomicNumber(atom);
        Ne = asfem_tools::getFixedElements(gridName,rN,order,atomicN);
        //std::cout<<"Points for this run: "<<Ne<<std::endl;
    }
    else{
        Ne = 100;
    }
    std::cout<<"Elements for this run: "<<Ne<<std::endl;
    //
    ASFEM asfem{Ne, order, poissNe, femModel,r0, rN, gridName,atomicN,atomicModel,integrals};
    asfem.printInitialData();
    asfem.startProgram();
    //asfem.printWfn(1);



    return 0;
}