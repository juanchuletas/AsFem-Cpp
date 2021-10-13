#if !defined(_FEM_FIXED_POINTS)
#define _FEM_FIXED_POINTS
#include "../Grids/Grid.hpp"
#include "../Tools/external_potential.hpp"
#include "../Tools/matrices_operations.hpp"
#include "../Tools/math_modules.hpp"
#include "../FEM_Matrices/fixed_points_analitic_vr_matrices.hpp"
#include "../FEM_Matrices/fixed_points_kinect_matrices.hpp"
#include "../FEM_Matrices/fixed_points_num_pot_matrices.hpp"
#include "../FEM_Matrices/fixed_points_overlap_matrices.hpp"
#include "../FEM_Matrices/fixed_points_mixed_vr_matrices.hpp"


class FixedPoints {

    //This class performs the fem fixed points model
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
        Grid<double> femgrid;
        double *lij{nullptr};
        double *bvec{nullptr};
        int poisson_domain_size;
        //Matrix<double> sij, kij, vij;//FEM MATRICES--
        //double *eSum{nullptr}, eSize{nullptr};
    private:
        void buildLinkMAtrix();

    public:

        FixedPoints();
        FixedPoints(int _Ne, int _order, Grid<double> &inputGrid);
        ~FixedPoints();
        void assambleMatrices(double *kij, double *sij, double *vij,double *boundaryVec,int atomicN);
        void femNumericalIntegration(double *inputMatrix,double *inputVector);
        void extractBCvector(double *l_mat, double *boundaryVec,int nodes);
        void reduceMatrix(double *matG, double *mat,int pts);
        int getBCSize();
        void setPoissonEquationData(double *matA, double *boundaryVector,int domain_size);
        double * getBoundaryConditionVector();
    protected:
        void solvePoissonEquation(double *xvector ,double *sourceVector);
        void getSourceVector(double *sourceVec,double *MatS, double *rightVector,double bcValue);
        


};



#endif // _FEM_FIXED_POINTS
