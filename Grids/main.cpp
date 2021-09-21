#include <iostream>
#include "Grid.hpp"
//#include "../Matrix/Matrix.hpp"


int main (){
    grid_tools::froese_fischer::hmf = 1.0/32.0;
    grid_tools::froese_fischer::rmf = 5.0;
    int Ne;
    int order = 2;
    double r0 = 0;
    double rc = 2.0;
    double rinfty = 250.0;
    int atomicN = 4.0;  
    double totpoints = grid_tools::froese_fischer::inverseKernel(rinfty,atomicN);
    int points = floor(totpoints);
    std::cout<<"Points gained: "<<points<<std::endl;
    if(points%2==0){
        printf("True\n");
        points = points - 1;
    }
    Ne = (points-1)/order;
    std::cout<<"Elements for this run: "<<Ne<<std::endl;
    std::cout<<"Points for this run: "<<Ne*order+1<<std::endl;
     std::string mesh = "Froese-Fischer";
    std::string mesh2 = "Chebyshev";
    std::string kind = "Fixed Elements";
    Grid<double> myGrid;
    myGrid.setGridData(r0,rinfty,Ne,order,mesh,atomicN);
    myGrid.createGrid(kind);
  
/*     int myIndex = myGrid.forcedInsertion(rc);
      for(int i=0; i<myGrid.size(); i++){
        printf("r[%d] = %lf\n",i, myGrid[i]);
    }
    myGrid.refineNear(myIndex,0.006,1); */
    for(int i=0; i<myGrid.size(); i++){
        printf("r[%d] = %lf\n",i, myGrid[i]);
    }
    /* int p = order+1;
    int k=0;
    for(int i=0; i<Ne; i++){
        for(int j=0; j<order+1; j++){
            printf("E[%d].n[%d].x = %d = %lf\n",i,j,k, myGrid[k]);
            k++; 
        }
        k--;
    }  */
    //myGrid.printGrid();


    return 0;
}