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
FEM<T>::FEM(int _Ne, int _order)
 :Ne{_Ne},order{_order}{

    globalSize = (Ne*order)+1;
    bcSize = globalSize-2;
}
template<class T>
FEM<T>::~FEM(){
    std::cout<<"FEM destructor called\n";
}
template<class T>
void FEM<T>::buildLinkMAtrix(){
	for(int i=0; i<Ne; i++)
        {
                for(int j=0; j<poly; j++)
                {
                        linkMat[poly*i+j] = poly*i+j-i;
                }
        }
}
template<class T>
void FEM<T>::setFemData(int inNe,int inOrder){
    order = inOrder;
    Ne = inNe;
    poly = order+1;
    linkmatSize = Ne*(order+1);
    linkMat.setMatrix(linkmatSize);
    buildLinkMAtrix();
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