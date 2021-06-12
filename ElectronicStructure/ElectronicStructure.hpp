#if !defined(_ELECTRONIC_STRUCTURE_H_)
#define _ELECTRONIC_STRUCTURE_H_
#include <iostream>
#include<string>

class ElectronicStructure : public FEMFP {
    
    private:
        int OccOrb;
        int numOfBasis;
        int charge;
        int atomicN;
    public:
        ElectronicStructure();
        void getExchangeDensity(double *rhox, double *wfn, int femBasis, int orb_a, int orb_b);
        void getDensity(double *rho, double *wfn, int femBasis, int occ);
        ~ElectronicStructure();
    /* namespace closedShell{

        void performHartreeFock();
    }

    namespace openShell{
        void performHartreeFock();
    } */


};

#endif // _ELECTRONIC_STRUCTURE_H_
