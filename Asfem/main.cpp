#include<iostream>
#include "Elements/searchelements.hpp"
#include "Asfem.cpp"

int main (){



    int Ne = 5;
    int order = 1;
    double r0 = 0.0;
    double rN = 20.0;
    std::string atom = "H";
    std::string atomicModel = "Free-Atom";
    std::string gridName = "Froese-Fischer";
    std::string femModel = "Fixed Elements";
    std::string integrals = "Numerical";
    Asfem<double> asfem{r0,rN,order,atom,gridName,femModel,atomicModel,integrals};
    asfem.printInputData();
    //asfem.runProgram();

    return 0;
}
