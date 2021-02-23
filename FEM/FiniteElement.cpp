#include "FiniteElement.hpp"



template<class T>
FEM<T>::FEM(){
    Ne = 0; 
    order = 0;
    globalSize = 0;
    poly = 0; 
    bcSize = 0; 
    linkmatSize = 0; 
}
template<class T>
FEM<T>::FEM(double _r0, double _rN, int _Ne, int _order,std::string _name)
 :Ne{_Ne}, order{_order},r0{_r0},rN{_rN} {

    globalSize = (Ne*order)+1;
    bcSize = globalSize-2;
    //femGrid.setGridData()

    /* eMatS = new T[elementMatSize];
    eMatK = new T[elementMatSize];
    eMatV = new T[elementMatSize]; */
    /* sij.setMatrix(bcSize);
    kij.setMatrix(bcSize);
    hij.setMatrix(bcSize);
    vij.setMatrix(bcSize);
    
    sij.fillWithZeroes();
    kij.fillWithZeroes();
    hij.fillWithZeroes();
    vij.fillWithZeroes(); */ 

}
template<class T>
void FEM<T>::setFemData(double _r0, double _rN, int _Ne, int _order,std::string meshName,int atomicN){
    order = _order;
    Ne = _Ne;
    poly = order+1;
    linkmatSize = Ne*(order+1);
    linkMat.setMatrix(linkmatSize);
    buildLinkMAtrix();
    femGrid.setGrid
    femGrid.createGrid(meshName,atomicN);
}
template<class T>
void FEM<T>::setFemData(int inNe,int inOrder,std::string meshName){
    order = inOrder;
    Ne = inNe;
    poly = order+1;
    linkmatSize = Ne*(order+1);
    linkMat.setMatrix(linkmatSize);
    buildLinkMAtrix();
    femGrid.createGrid(meshName);
}
template<class T>
int & FEM<T>::getLinkMatIndex(int i){
 return linkMat[i];
}
template<class T>
T *FEM<T>::getOverlap(){
    return ElementalOverlapMatrix(order);
}
template<class T>
T *FEM<T>::getKinect(){
    return ElementalKinectMatrix(order);
}
template<class T>
T *FEM<T>::getPotential(double *vij){
    return ElementalPotentialMatrix(order,vij);
}
// *** END METHODS