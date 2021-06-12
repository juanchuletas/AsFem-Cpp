#include "ElectronicStructure.hpp"



ElectronicStructure::ElectronicStructure(){
    
}
void ElectronicStructure::getExchangeDensity(double *rhox, double *wfn, int femBasis, int orb_a, int orb_b){
    for(int i=0; i<femBasis; i++){
        rhox[i] = 0.0;
        rhox[i] = wfn[orb_a*femBasis + i]*wfn[orb_b*femBasis + i];
    }
    



    /* void getDensity(double *rho, double *wfn, int fembasis,int occ=2.0){
    }
    namespace closedShell{


    }
    namespace openShell{


    } */




