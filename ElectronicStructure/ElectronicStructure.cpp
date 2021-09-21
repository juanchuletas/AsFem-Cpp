#include "ElectronicStructure.hpp"


namespace electronicStructure{

    //Constructors: ******************
    closedShell::closedShell()
    :nucleii{0},charge{0},bcnodes{0},
        gridnodes{0},orbitals{0},numElectrons{0}{

    }
    closedShell::closedShell(int _nuc,int _char)
    :nucleii{_nuc},charge{_char}
    {
        std::cout<<"Closed-Shell Constructor\n";
        numElectrons = nucleii-charge;
        orbitals  = numElectrons/2;
    }
    // DESTRUCTOR:
    closedShell::~closedShell(){
    }
    //******************************************
     //Closed Shell Methods:***********************
    void closedShell::setNumOfGridPoints(int _gridPoints){
        gridnodes = _gridPoints;
        bcnodes = _gridPoints-2;
    }
    void closedShell::computeDensityMatrix(double *densMat,double *wfn){
        int k=0; 
        for(int i=0; i<bcnodes; i++){
            for(int j=0; j<bcnodes; j++){
                densMat[k] = 0.0;
                for(int orb=0; orb<orbitals; orb++){
                    densMat[k] += 2.0*(wfn[i + orb*bcnodes]*wfn[j + orb*bcnodes]);
                }
                //densMat[k] = 2.0*densMat[i];
                /* if(i==j){
                    printf("rho[%d] = %lf\n", k,densMat[k]);
                } */
                k++; 

            }
        }
    }
    double closedShell::solveHartreeFock(double *fij, double *hij, double *densmat){
        double energy;
        for(int i=0; i<bcnodes*bcnodes; i++){
            energy = energy + 0.5*(densmat[i]*hij[i]);
            energy = energy + 0.5*(densmat[i]*fij[i]);
        }        
        return 0.5*energy;
    }




}



