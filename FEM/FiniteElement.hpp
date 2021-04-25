#if !defined(_FEM_H_)
#define _FEM_H_
#include<iostream>
#include<string>
#include "../Grids/Grid.hpp"
#include "../Tools/math_modules.hpp"
#include "../FEM_Matrices/fixed_elements_num_potential_matrices.hpp"
#include "../FEM_Matrices/fixed_elements_overlap_matrices.hpp"
#include "../FEM_Matrices/fixed_elemets_kinect_matrices.hpp"
#include "../FEM_Matrices/fixed_points_analitic_vr_matrices.hpp"
#include "../FEM_Matrices/fixed_points_kinect_matrices.hpp"
#include "../FEM_Matrices/fixed_points_num_pot_matrices.hpp"
#include "../FEM_Matrices/fixed_points_overlap_matrices.hpp"

class FEM {
    
    protected: 
        int Ne, order;
        int globalSize;
        int poly;
        int bcSize;
        int linkmatSize;
        std::string femModel;
        std::string gridType;
        int poissNe;
        int poissNodes;
        //double r0,rN;
        Matrix<int> linkMat;
        Matrix<double> elsize;
        Grid<double> femgrid;
        Matrix<double> sij, kij, vij;//FEM MATRICES--
        //double *eSum{nullptr}, eSize{nullptr};
        double *bvec{nullptr};
    private:
        void buildLinkMAtrix();
        void applyBoundaryConditions();
        void reduceMatrix(double *, Matrix<double> &mat, int size);
        void reduceMatrix(double *matG, double *mat,int size);
        void extractBCvector(double *l_mat, int nodes);

        //T* OverlapElementalMatrices();

    public:
        FEM(); 
        FEM(int _Ne, int _order, std::string _femModel, std::string _gridName); 
        FEM(int _Ne, int _order,int _poissonNe, std::string _femModel,std::string gridName);
        ~FEM();
        // **** METHODS**********
        void solvePoissonEquation(double *hpot,double *rho_r,double hp); //Working! Don't touch
        void assamblePoissonMatrices(double *lij, double *uij,int pNe);
        void buildFemGrid(int atomicN,double r0,double rN); //Working
        void setFemData(int inNe, int inOrder);
        int& getLinkMatIndex(int i);
        // ******* FIXED ELEMENTS MODEL METHODS***********

        void assambleMatricesFixedElements(double *v); // for Kinect, Overlap adnd Potential Matrices
        void fixedElementsNumIntegration(double *mat,double *vec);
        //****** FIXED POINTS METHODS  ****************
        void assambleMatricesFixedPoints(int atomicN); 
        void fixedPointsNumIntegration(double *mat, double *vec);//Full numeric integration with no equal spaced elements

};


#include "FiniteElement.cpp"

#endif // _FEM_H_

