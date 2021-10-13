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
template<class T>
Grid<T>& Grid<T>::operator=(const Grid<T> &source){
       grid = source.grid;
}
// ***** END OPERATORS************
// *************** METHODS ************************
template<class T>
void Grid<T>::setSize(int _size){
    totalNodes = _size;
    grid.setMatrix(_size);
}
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
void Grid<T>::createGrid(std::string name)
{
    if(meshType=="Froese-Fischer"){
        buildAtomic(atomicN,name);
    }
    else if(meshType=="Chebyshev"){
        buildChebyshev();
    }
}
template<class T>
int Grid<T>::forcedInsertion(double cutRad){
    //Only working for Froese-Fischer Grid
    std::cout<<"Force Insertion Module\n";
    /* double target;
    int index;
    for(int i=0; i<this->size(); i++){
        index = i;
        target = grid[i];
        if(cutRad<target){
            if(index%2!=0){
                index++;
                grid[index] = cutRad;
                target = grid[index];
            }
            else{
                grid[i] = cutRad;
            }
            break;
        }
    } */
    //printf("DELETED ITEM: %lf IN POSITION %d\n",target,index);
    double rcPoints  = grid_tools::froese_fischer::inverseKernel(cutRad,atomicN);
    int index = floor(rcPoints);
    std::cout<<"inverse value: "<<rcPoints<<std::endl;
    //nt points = index;
    switch (order)
    {
    case 2:
            if(index%2==0){
                //index = index+2;
                grid[index] = cutRad;
            }
            else{
                index = index+1;
                grid[index] = cutRad;
            }
        break;
    case 3:
            if(index%order!=0){
                if(index%2==0){
                    index++;
                }
                else{
                    index--;
                }
                grid[index] = cutRad;
            }
            else{
                grid[index] = cutRad;
            }
        break;
    }
    int points = index;
    printf("RC index = %d\n",points);
    printf("Elements at RC = %d\n",(points-1)/order);
    int poiss_tot_nodes = points+1;
    int poiss_bc = poiss_tot_nodes-2;
    std::cout<<"Poisson equation total points: "<<index+1<<std::endl;
    std::cout<<"Poisson equation with BC points: "<<poiss_bc<<std::endl;


    return points;
}
template<class T>
void Grid<T>::refineNear(int point, double delta, int steps){
    int index = point;
    int target = (index-1)/order;
    std::cout<<"Refining: "<<steps<<" steps before and after: "<<index<<std::endl;
    for(int i=-steps; i<=steps; i++){
        int j=index+i;
        //double h = grid[j]-grid[j-1];
       // double new_h = h/5.0;
        printf("grid[%d] = %lf\n",j,grid[j]);
        grid[j] = grid[index] + delta*i;
        printf("grid[%d] = %lf\n",j,grid[j]);
        /* printf("grid[%d] = %lf + %lf\n",j,grid[j],new_h);
        if(i==0){
            grid[j] = grid[j] + new_h*0.0;
        }
        else{
            grid[j] = grid[j] - delta*i;
        }
        printf("grid[%d] = %lf\n",j,grid[j]) */;
    }
}
template<class T>
void Grid<T>::chebyshev(){
    double xi,ra,rb;
    double  rMax;
    double xmax = 1. - cos( totalNodes*M_PI / (2.*(totalNodes + 1.)));
    T *xprev;
    xprev = new T[Ne];
    T faca = (rN - r0)/xmax;
    T facb = (r0*xmax)/xmax;
    for ( int i = 0; i < totalNodes; i++){
        xi = 1. - cos(i*M_PI / (2.*( totalNodes + 1.)));
        grid[i] = faca * xi + facb;
        //printf("%lf\n",xprev[i]);
    }
    grid[totalNodes-1] = rN;
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
    grid[0] = 0.f;
    grid[totalNodes-1] = rN;
    //printf("count = %d\n",count);

    delete [] rinitial;

}
template<class T>
void Grid<T>::buildAtomic(int atomicN, std::string name){
    int totnodes = Ne*order+1;
    
    double ri,f_r1,f_r2;
    double nucleii = static_cast<double>(atomicN);
   

    double rmin  = grid_tools::froese_fischer::kernel(0,nucleii);
    double rmax  = grid_tools::froese_fischer::kernel(totnodes,nucleii);
    
    f_r1 = (rN-r0)/(rmax-rmin);
    f_r2 = (r0*rmax-rN*rmin)/(rmax-rmin);   
    for(int i=1;i<totnodes;i++)
    {
        //ri = FroeseFischer(i,atomicN,rmf,hmf);
        ri = grid_tools::froese_fischer::kernel(i,nucleii);
        //grid[i] = f_r1*ri + f_r2; 
        grid[i] = ri;
        //printf("Vertex value x[%d] = %lf\n",i,x[i]);
    }
    grid[0] = 0.f;
    grid[totnodes-1] = rN;
}
template<class T>
void Grid<T>::froeseFischer(int atomicN){
    int totnodes = Ne*order+1;
    
    double ri,f_r1,f_r2;
    double nucleii = static_cast<double>(atomicN);
   

    double rmin  = grid_tools::froese_fischer::kernel(0,nucleii);
    double rmax  = grid_tools::froese_fischer::kernel(totnodes,nucleii);
    
    f_r1 = (rN-r0)/(rmax-rmin);
    f_r2 = (r0*rmax-rN*rmin)/(rmax-rmin);   
    for(int i=1;i<totnodes;i++)
    {
        //ri = FroeseFischer(i,atomicN,rmf,hmf);
        ri = grid_tools::froese_fischer::kernel(i,nucleii);
        //grid[i] = f_r1*ri + f_r2; 
        grid[i] = ri;
        //printf("Vertex value x[%d] = %lf\n",i,x[i]);
    }
    grid[0] = 0.f;
    grid[totnodes-1] = rN;
}
template<class T>
void Grid<T>::exponential(){
    double denom = pow(totalNodes,2);
    double nume;
    double fracc;
    double basis = (1+rN);
    for(int i=1; i<totalNodes; i++){
        nume = pow(i,2);
        fracc = nume/denom;
        grid[i] = pow(basis,fracc) - 1;
    }
    grid[0] = 0.0;
    grid[totalNodes-1] = rN;
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
    //grid.printMatrix();
    for(int i=0; i<totalNodes; i++){
        printf("r[%d] = %lf\n",i,grid[i]);
    }
}
// *************** END METHODS ********************
