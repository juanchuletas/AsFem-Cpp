#if !defined(_ASFPM_H_)
#define _ASFPM_H_
#include <fstream>
#include "../FEM_FP/FiniteElementFixedPoints.hpp" //Only add the 
#include "../Tools/external_potential.hpp"
#include "../Tools/matrices_operations.hpp"
#include "../Tools/math_modules.hpp"
#include "../ElectronicStructure/ElectronicStructure.hpp"
    
class ASFPM : public FEMFP, public electronicStructure::closedShell {

    //This class performs the Finite Elemennt Method for Atomic Structure under the 
    // Fixed Points Model
    // CLASS MEMBERS
    
        int occOrb;
        double *wfn{nullptr}, *eigenVal{nullptr};
        double r0;
        double rN;
        int atomicN,angular;
        double cutRad, wallValue, lambda;
        int total_nodes;
        int fem_nodes;
        int charge;
        int numElectrons; 
        double totQ;
        std::string atomName;
        std::string atomicModel;
        std::string confType;
        std::string integrationScheme;
        double *rho{nullptr},*vr{nullptr}; 
    // VERY IMPORTANT MATRICES:
    //Matrix<double> sij, vij, kij;
    //double *uij{nullptr}, *lij{nullptr}; //Only when the poisson problem has a different size 
    //Private Methods
        void getDensity(const double ocupation=2.0);
        
        void getExternalPotential(); //get the -z/r
        void performSCF();
        double* computeHartreePotential(double *hpot,int rcIndex);
        void getExternalPotentialMatrix();
    
   
    double computeTotalCharge(int rcindex);
    //double integrateElement();

    public:
        ASFPM();
        ASFPM(int _Ne, int _order,std::string _grid,int _atomicN,int _charge, int _angular,double _rN);
        ~ASFPM();
        //void setData(std::string _atomicmodel);
        void getExchangeDens(double *rox, int a, int b);
        void setData(std::string _atomicmodel, double _rc=0.0, double _wall=0.0);
        void startProgram(); //Runs the program
        void load();
        void singleDiagonalization(double *hij);
        double *getOrbitals(int numOfOrb);
        double *exchangeAuxPotential();
        double *hartreePotential();
        void getDensity(double *rhor);
        void sfc(double *fij, double *hij, double *vhij, double *vxij,double &energy0);
        double* computeHartreePotential();
        double* computeExchangePotential();
        void getDensityMatrix(double *densMat);
        double energyHF(double *hij, double *fij,double *densMat);
        void wfnNormalization();
        void printWfn(int orb );
        double computeBoundCond();
};
double integrateElement(int ei,int order,double *feMatS, int *link_mat, double coeff, double *cf);

#endif // _ASFPM_H_