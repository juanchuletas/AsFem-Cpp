#if !defined(_FEM_H_)
#define _FEM_H_
#include<iostream>
#include<string>
#include "../Grids/Grid.hpp"
#include "overlap_elemental_matrices.hpp"
#include "kinect_elemental_matrices.hpp"
#include "potential_elemental_matrices.hpp"


template<class T> class FEM {

    int Ne, order;
    int globalSize;
    int poly;
    int bcSize;
    int linkmatSize;
    //double r0,rN;
    Matrix<int> linkMat;
    private:
        void buildLinkMAtrix();
        //T* OverlapElementalMatrices();

    public:
        FEM(); 
        FEM(int _Ne, int _order);
        ~FEM();
        void setFemData(int inNe, int inOrder);
        int& getLinkMatIndex(int i);
        //void buildElementalMatrices();
        T* getOverlap();
        T* getKinect();
        T* getPotential(double *vij);


};

#include "FiniteElement.cpp"


#endif // _FEM_H_

