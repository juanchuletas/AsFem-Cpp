#if !defined(_ASFEM_H_)
#define _ASFEM_H_
#include "../FEM/FiniteElement.hpp"
#include "Tools/searchelements.hpp"
#include "Tools/external_potential.hpp"
#include "Tools/math_tools.hpp"
#include "Tools/matrices_operations.hpp"
#include "Tools/math_modules.hpp"
#include <fstream>

class ASFEM: public FEM {

    // CLASS MEMBERS
    double r0;
    double rN;
    int atomicN,angular;
    double cutRad, wallValue, lambda;
    int total_nodes;
    int fem_nodes;
    int occOrb;
    int charge; 
    std::string atomName;
    std::string atomicModel;
    std::string confType;
    std::string integrationScheme;
    double *vr{nullptr}, *vh{nullptr};
    double *wfn{nullptr}, *eigenVal{nullptr},*rho{nullptr}; 
    // VERY IMPORTANT MATRICES:
    //Matrix<double> sij, vij, kij;
    double *uij{nullptr}, *lij{nullptr}; //Only when the poisson problem has a different size 
    //Private Methods
    void computeHartreePotential();
    void getDensity();
    void diagonalize();
    void getExternalPotential(); //get the -z/r 
    void HartreeFock();
    void performSCF();
    void getExternalPotentialMatrix();
    void singleDiagonalization();
    void wfnNormalization();
    void wfnNormalization(std::string name);
    //double integrateElement();


    public:
        ASFEM();
        ASFEM(int _Ne, int _order, int _poissonNe, std::string _femModel,double rin, double rInfty,std::string gridType,int atom,std::string intAtomicModel,std::string integrationType);
        ~ASFEM();
        //METHODS
        void startProgram();
        void printInitialData();
        void printWfn(int orbital);
        //Quantun Chemistry Methods
        



}; //_ASFEM_H_
double integrateElement(int ei,int order,double *feMatS, int *link_mat, double coeff, double *cf){
    int poly = order +1; 
    double phi[poly];
    for(int j=0; j<poly; j++){
        phi[j] = cf[link_mat[poly*ei+j]];
    }
    double value = 0.0;
    for(int mu=0; mu<poly; mu++){
        for(int nu=0; nu<poly; nu++){
            value = value + phi[mu]*phi[nu]*(coeff*feMatS[poly*mu + nu]);
        }
    }
    return value;

}
#endif