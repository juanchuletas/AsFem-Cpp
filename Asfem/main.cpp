#include<iostream>
#include "Asfem.cpp"

int main (){
    int Ne = 5;
    int order = 1;
    double r0 = 0.0;
    double rN = 10.0;
    int atom = 2;
    std::string name = "Froese-Fischer";
    Asfem<double> femito{r0,rN,Ne,order,atom,name};
    femito.runProgram();

    return 0;
}