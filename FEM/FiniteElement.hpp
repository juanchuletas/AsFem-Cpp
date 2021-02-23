#if !defined(_FEM_H_)
#define _FEM_H_
#include<iostream>
#include<string>
#include "../Grids/Grid.cpp"
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

    protected: 
        Grid<T> femGrid;
    private:
        void buildLinkMAtrix();
        //T* OverlapElementalMatrices();

    public:
        FEM(); 
        FEM(double _r0, double _rN, int _Ne, int _order,std::string _name);
        ~FEM();
        void setFemData(double _r0, double _rN, int _Ne, int _order,std::string _name,int atomicN);
        void setFemData(int inNe, int inOrder,std::string name);
        int& getLinkMatIndex(int i);
        //void buildElementalMatrices();
        T* getOverlap();
        T* getKinect();
        T* getPotential(double *vij);


};




#endif // _FEM_H_

