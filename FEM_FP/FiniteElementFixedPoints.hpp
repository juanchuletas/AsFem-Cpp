#if !defined(_FEM_FIXED_POINTS)
#define _FEM_FIXED_POINTS  
#include <iostream>
#include "../Grids/Grid.hpp"
#include "../Tools/math_modules.hpp"
#include "../FEM_Matrices/fixed_points_analitic_vr_matrices.hpp"
#include "../FEM_Matrices/fixed_points_kinect_matrices.hpp"
#include "../FEM_Matrices/fixed_points_num_pot_matrices.hpp"
#include "../FEM_Matrices/fixed_points_overlap_matrices.hpp"
#include "../FEM_Matrices/fixed_points_mixed_vr_matrices.hpp"

class FEMFP {

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
        void divideOverGridPoints(double *inVector);
        void reduceMatrix(double *, Matrix<double> &mat, int size);
        void reduceMatrix(double *matG, double *mat,int size);
        void reduceMatrix(double *matG, double *mat,int sized,int pts);
        void extractBCvector(double *l_mat, int nodes);
        void extractBCvector(double *l_mat,double *bcVec,int nodes);
    public:
        FEMFP(); 
        FEMFP(int _Ne, int _order, std::string _femModel, std::string _gridName);
        ~FEMFP();
        void assambleMatricesFixedPoints(int atomicN); 
        void assambleMatricesFixedPoints(int atomicN,int points); //Only for special ocations 
        void fixedPointsNumIntegration(double *mat, double *vec);//Full numeric integration with no equal spaced elements
        void buildFemGrid(int atomicN,double r0,double rN); //Working
        // POISSON SOLVER
        void solvePoissonEquation(double *hpot,double *rho_r,double hp); //Working! Don't touch
        void assamblePoissonMatrices(double *lij, double *uij,double *bcVec,int rcIndex);
};



#endif // _FEM_FIXED_POINTS
