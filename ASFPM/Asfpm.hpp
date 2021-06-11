#if !defined(_ASFPM_H_)
#define _ASFPM_H_
#include "../FEM_FP/FiniteElementFixedPoints.hpp" //Only add the 
#include "../Tools/external_potential.hpp"
#include "../Tools/matrices_operations.hpp"
#include "../Tools/math_modules.hpp"
#include "../ElectronicStructure/ElectronicStructure.hpp"
    
class ASFPM : public FEMFP{

    //This class performs the Finite Elemennt Method for Atomic Structure under the 
    // Fixed Points Model
    // CLASS MEMBERS
    double r0;
    double rN;
    int atomicN,angular;
    double cutRad, wallValue, lambda;
    int total_nodes;
    int fem_nodes;
    int occOrb;
    int charge;
    int numElectrons; 
    double totQ;
    std::string atomName;
    std::string atomicModel;
    std::string confType;
    std::string integrationScheme;
    double *wfn{nullptr}, *eigenVal{nullptr},*rho{nullptr},*vr{nullptr}; 
    // VERY IMPORTANT MATRICES:
    //Matrix<double> sij, vij, kij;
    //double *uij{nullptr}, *lij{nullptr}; //Only when the poisson problem has a different size 
    //Private Methods
    void getDensity(const double ocupation);
    void getExternalPotential(); //get the -z/r
    void performSCF();
    double* computeHartreePotential(double *hpot);
    double* computeExchangePotential(double *xpot);
    double* computeHartreePotential(double *hpot,int rcIndex);
    void getExternalPotentialMatrix();
    void wfnNormalization();
    void getDensityMatrix(double *densMat);
    double energyHF(double *hij, double *fij,double *densMat);
    double computeTotalCharge(int rcindex);
    //double integrateElement();

    public:
        ASFPM();
        ASFPM(int _Ne, int _order,std::string _grid,int _atomicN,int _charge, int _angular,double _rN);
        ~ASFPM();
        void setData(std::string _atomicmodel);
        void setData(std::string _atomicmodel, double _rc, double _wall);
        void startProgram(); //Runs the program
        void load();
        void singleDiagonalization();
        double *getOrbitals(int numOfOrb);

};
double integrateElement(int ei,int order,double *feMatS, int *link_mat, double coeff, double *cf);

#endif // _ASFPM_H_
