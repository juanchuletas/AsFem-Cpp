#if !defined(_ATOMIC_STRUCTURE)
#define _ATOMIC_STRUCTURE
#include <iostream>
#include <fstream>
#include "../FP_Model/FemFixedPoints.hpp"
class Atomic : public FixedPoints{
    
    protected:
        int numOrb;
        int virtualOrbs;
        //double *wfn{nullptr}, *eigenVal{nullptr};
        double r0;
        double rN;
        int atomicN,angular;
        double cutRad, wallValue, lambda;
        int total_nodes;
        int fem_nodes;
        int numElec; 
        double totQ;
        std::string atomName;
        std::string atomicModel;
        std::string confType;
        std::string integrationScheme;
        double *orbVec;
        //double *rho{nullptr}; 
    public:
        // ---- Constructors and destructors
        Atomic();
        Atomic(int _Ne, int _order,Grid<double> _grid, int atomicN, int _numElec,int _orbs ,int virtual_orbs,double _rMax);
        ~Atomic();
        void printOrbital(std::string name,int orbOfinterest);
        void printWfn();
        // --- Methods
    protected:
        void wfnNormalization(double *wfn);
        int orbitalPhase(int orbOfinterest);
        double * divideOverGridPoints(double *input);
        void integrateHartreePotential();
        void integrateExchangePotential();
        void getOrbitals(double *matCoeff);
        double energyHF(double *hcore, double *fock_mat,double *densmat);
        void rayleighQuotient(double *hcoremat, double *smat,double *exchagevec);
        void samePhase(int orb1, int orb2);

};
double integrateElement(int ei,int order,double *feMatS, int *link_mat, double coeff, double *cf);
#endif // _ATOMIC_STRUCTURE
