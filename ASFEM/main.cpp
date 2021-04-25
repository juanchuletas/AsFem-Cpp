#include<iostream>
#include<string>
#include "Asfem.cpp"



int main ()
{

    grid_tools::froese_fischer::rmf = 5.0;
    grid_tools::froese_fischer::hmf = 1.0/24;
    int Ne;
    int order = 2;
    double r0 = 0.0;
    double rN = 20.0;
    int poissNe = 100;
    int angular = 0;
    int charge = 0;
    double lambda = 0.f, wallValue = 0.f, rC = 0.f; 
    std::string atom = "He";
    std::string atomicModel = "Free-Atom";
    std::string gridName = "Froese-Fischer";
    std::string femModel = "Fixed Points";
    std::string integrals = "Analitic";
    std::string confType = "Free";
    int atomicN = getAtomicNumber(atom);
    if(femModel=="Fixed Points"){
        printf("Fixed Points\n");
        int atomicN = getAtomicNumber(atom);
        Ne = asfem_tools::getFixedElements(gridName,rN,order,atomicN);
        //std::cout<<"Points for this run: "<<Ne<<std::endl;
    }
    else{
        Ne = 200;
    }
    std::cout<<"Elements for this run: "<<Ne<<std::endl;
    /*ASFEM ASFEM(std::string _femModel, int Ne, 
    int order,std::string _atomicModel, 
    double _lambda,std::string _confType,
    double _Rc,double _wallVal,int _atomicN, 
    int _charge,int _angular, double rInfty,
    std::string _integrals);*/
    ASFEM asfem{femModel,Ne,order,atomicModel,lambda,confType,rC,wallValue,atomicN,charge,angular,gridName,rN, integrals};
    asfem.startProgram();
    asfem.printInitialData();
    /* ASFEM asfem{Ne, order, poissNe, femModel,r0, rN, gridName,atomicN,atomicModel,integrals};
    asfem.printInitialData();
    asfem.startProgram(); */
    asfem.printWfn(1);



    return 0;
}