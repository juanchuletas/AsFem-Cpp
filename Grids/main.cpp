#include <iostream>
#include "Grid.hpp"
//#include "../Matrix/Matrix.hpp"


int main (){
    grid_tools::froese_fischer::hmf = 1.0/24.0;
    grid_tools::froese_fischer::rmf = 5.0;
    int Ne = 191;
    int order = 2;
    double r0 = 0;
    double rc = 1.5;
    double rinfty = 20.0;
    int atomicN = 2;  
    std::string mesh = "Froese-Fischer";
    std::string mesh2 = "Chebyshev";
    std::string kind = "Fixed Elements";
    Grid<double> myGrid;
    myGrid.setGridData(r0,rinfty,Ne,order,mesh,2);
    myGrid.createGrid(kind);
    myGrid.forceInsertion(rc);
    printf("Grid size: %d\n",myGrid.size());
    for(int i=0; i<myGrid.size(); i++){
        printf("r[%d] = %lf\n",i, myGrid[i]);
    }
   
    //myGrid.printGrid();


    return 0;
}