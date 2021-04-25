#if !defined(_ASFPM_H_)
#define _ASFPM_H_
#include "../FEM/FiniteElement.hpp"
class ASFPM : public FEM {

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
    double *vr{nullptr}, *vh{nullptr};
    double *wfn{nullptr}, *eigenVal{nullptr},*rho{nullptr}; 
    // VERY IMPORTANT MATRICES:
    //Matrix<double> sij, vij, kij;
    double *uij{nullptr}, *lij{nullptr}; //Only when the poisson problem has a different size 
    //Private Methods
    void getDensity();
    void getExternalPotential(); //get the -z/r
    void performSCF();
    double* computeHartreePotential(double *hpot);
    void getExternalPotentialMatrix();
    void singleDiagonalization();
    void wfnNormalization();
    void getDensityMatrix(double *densMat);
    double energyHF(double *hij, double *fij,double *densMat);
    //double integrateElement();

    public:
        ASFPM();
        ASFPM(std::string femModel, int Ne, int order,std::string _atomicModel, double _lambda,std::string _confType,double _Rc,double _wallVal,int _atomicN, int _charge,int _angular, std::string gridType,double rInfty,std::string _integrals);
        
};

#endif // _ASFPM_H_
