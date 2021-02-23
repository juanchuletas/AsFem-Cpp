#include "Grid.hpp"




//********* CONSTRUCTORS AND DESTRUCTORS ***********
template<class T>
Grid<T>::Grid()
:Ne{0},order{0},totalNodes{0},polynomial{0},bcNodes{0},r0{0},rN{0}{
    meshType  = "NULL";
    atomicN = 0;
    //PENDEJOOOOOO

}
template<class T>
Grid<T>::Grid(double rinit,double rfinal,int inNe, int inOrder,std::string inType)
:Ne{inNe},order{inOrder},meshType{inType},r0{rinit},rN{rfinal}{
    totalNodes = Ne*order+1;
    bcNodes= totalNodes-2;
    polynomial = order+1;
    atomicN = 0;
    grid.setMatrix(totalNodes);
    elSize.setMatrix(Ne);
    printf("Grid size: %d\n",grid.getSize());
    printf("Elment size: %d\n",elSize.getSize());

}
template<class T>
Grid<T>::Grid(double rinit,double rfinal,int inNe, int inOrder,std::string inType,int inatomicN)
:Ne{inNe},order{inOrder},meshType{inType},r0{rinit},rN{rfinal},atomicN{inatomicN}{
    totalNodes = Ne*order+1;
    bcNodes= totalNodes-2;
    polynomial = order+1;
    grid.setMatrix(totalNodes);
    elSize.setMatrix(Ne);
    printf("Grid size: %d\n",grid.getSize());
    printf("Elment size: %d\n",elSize.getSize());

}


template<class T>
Grid<T>::~Grid(){

}
//********* END CONSTRUCTORS AND DESTRUCTORS ***********
// **** OPERATORS ************
template<class T>
T& Grid<T>::operator[](int i){
    return grid[i];
}

// ***** END OPERATORS************
// *************** METHODS ************************
template<class T>
void Grid<T>::createGrid()
{
    if(meshType=="Froese-Fischer"){
        buildAtomic(atomicN);
    }
    else if(meshType=="Chebyshev"){
        buildChebyshev();
    }
}
template<class T>
void Grid<T>::buildChebyshev(){
    double xi,ra,rb;
    double  rMax;
    double xmax = 1. - cos( Ne*M_PI / (2.*(Ne + 1.)));
    T *xprev;
    xprev = new T[Ne];
    T faca = (rN - r0)/xmax;
    T facb = (r0*xmax)/xmax;

    for ( int i = 0; i < Ne; i++){
        xi = 1. - cos(i*M_PI / (2.*( Ne + 1.)));
        xprev[i] = faca * xi + facb;
        //printf("%lf\n",xprev[i]);
    }

  for(int i=0; i<Ne-1;i++){
        elSize[i] = xprev[i+1]-xprev[i];
    }
    elSize[Ne-1] = rN-xprev[Ne-1];
    int count =0;
    for(int i=0; i<Ne; i++){
        grid[i*order] = xprev[i];
        for(int j=0; j<order; j++){
            grid[i*order+j] = grid[i*order] + (j*elSize[i])/(order);
            count++;
        }
    }
    grid[totalNodes-1] = rN; 
    delete [] xprev;
}
template<class T>
void Grid<T>::buildAtomic(int atomicN){

    int totnodes = Ne*order+1;
    double rmf = 5.0;
    double hmf = 1.0/32.0;
    double ri,f_r1,f_r2;
    T *rinitial;
    rinitial = new T[Ne];

    double rmin = FroeseFischer(0,atomicN,rmf,hmf);
    double rmax = FroeseFischer(Ne,atomicN,rmf,hmf);
    
    f_r1 = (rN-r0)/(rmax-rmin);
    f_r2 = (r0*rmax-rN*rmin)/(rmax-rmin);   
    for(int i=0;i<Ne;i++)
    {
        ri = FroeseFischer(i,atomicN,rmf,hmf);
        rinitial[i] = f_r1*ri + f_r2;
        //printf("Vertex value x[%d] = %lf\n",i,x[i]);
    }

    for(int i=0; i<Ne-1;i++){
        elSize[i] = rinitial[i+1]-rinitial[i];
    }
    elSize[Ne-1] = rN-rinitial[Ne-1];
    int count =0;
    for(int i=0; i<Ne; i++){
        grid[i*order] = rinitial[i];
        for(int j=0; j<order; j++){
            grid[i*order+j] = grid[i*order] + (j*elSize[i])/(order);
            count++;
        }
    }
    grid[totalNodes-1] = rN;
    printf("count = %d\n",count);

    delete [] rinitial;

}
template<class T>
void Grid<T>::setGridData(double rinit, double rfinal,int inNe, int inorder,std::string inmeshType){
    r0 = rinit; rN = rfinal; Ne = inNe;
    order = inorder; meshType = inmeshType;
    totalNodes = Ne*order+1;
    bcNodes= totalNodes-2;
    polynomial = order+1;

    grid.setMatrix(totalNodes);
    elSize.setMatrix(Ne);
}
template<class T>
void Grid<T>::setGridData(double rinit, double rfinal,int inNe, int inorder,std::string inmeshType,int inatomicN){
    r0 = rinit; rN = rfinal; Ne = inNe;
    order = inorder; meshType = inmeshType;
    atomicN = inatomicN;
    totalNodes = Ne*order+1;
    bcNodes= totalNodes-2;
    polynomial = order+1;
    grid.setMatrix(totalNodes);
    elSize.setMatrix(Ne);
}
template<class T>
int Grid<T>::size(){
    return totalNodes;
}
template<class T>
T& Grid<T>::getGridItem(int i){

    return grid[i];
}
template<class T>
T& Grid<T>::getElementSize(int i){
    
    return elSize[i]; 
}
template<class T>
void Grid<T>::printGrid(){
    grid.printMatrix();
}
// *************** END METHODS ********************
