#include <iostream>
#include "Grid.hpp"
#include "../Matrix/Matrix.hpp"


int main (){

    int Ne = 10;
    int order = 3;
    double r0 = 0;
    double rc = 4.0;
    double rinfty = 10.0;
    int atomicN = 2;  
    std::string mesh = "Froese-Fischer";
    std::string mesh2 = "Chebyshev";
    Grid<double> myGrid,mygrid2;
    mygrid2.setGridData(rc,rinfty,Ne,order,mesh2,2);
    myGrid.setGridData(r0,rc,Ne,order,mesh,2);
    mygrid2.createGrid();
    myGrid.createGrid();
    for(int i=0; i<myGrid.size(); i++){
        printf("%lf\n", myGrid[i]);
    }
    for(int i=0; i<myGrid.size(); i++){
        printf("%lf\n", mygrid2[i]);
    }
    //myGrid.printGrid();


    return 0;
}