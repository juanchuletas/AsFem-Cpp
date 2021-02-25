#include <iostream>
#include "Grid.hpp"
#include "../Matrix/Matrix.hpp"


int main (){

    int Ne = 10;
    int order = 1;
    double r0 = 0;
    double rN = 10;
    int atomicN = 2;  
    std::string mesh = "Froese-Fischer";
    Grid<double> myGrid;
    myGrid.setGridData(r0,rN,Ne,order,mesh,2);
    myGrid.createGrid();
    for(int i=0; i<myGrid.size(); i++){
        printf("%lf\n", myGrid[i]);
    }
    //myGrid.printGrid();


    return 0;
}