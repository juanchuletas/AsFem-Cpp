#include<iostream>
#include<cmath>
#include "../ASFPM/Asfpm.hpp"
#include "../UserInput/UserInputData.hpp"
#include "../Tools/searchelements.hpp"
#define MIN_REQUIRED 2



int setNumberOfElements(double rN, int atomicN,int order){
    double totpoints = grid_tools::froese_fischer::inverseKernel(rN,atomicN);
    int points = floor(totpoints);
    if(points%2==0){
        points = points - 1;
    } 
    int Ne = (points-1)/order;

    return Ne;
}
int help() {
   printf("Usage: asfem.exe [-s ] [-n ] [-true]\n");
   printf("\t-s: a string a\n");
   printf("\t-n: a number\n");
   printf("\t-true: a single parameter\n");

   return 1;
}
//Electronic Structure namespace:
using namespace electronicStructure;
//*********************************


int main(int argc, char **argv){

    std::string atomSym,meshType,atomicModel,femModel,confType, integrals,multi;
    int angular,Ne,order,charge;
    double r0, rN,lambda, rC,wallValue;
    grid_tools::froese_fischer::rmf = 5.0;
    grid_tools::froese_fischer::hmf = 1.0/32.0;
    if(argc<MIN_REQUIRED){
        std::cout<<"IDIOT, YOU NEED AN INPUT FILE"<<std::endl;
        return help();
    }
    else{
        int flag = read_data(argv,femModel,Ne, order, atomicModel,confType,rC,wallValue,lambda,charge,atomSym,multi,angular,r0,rN,meshType,integrals);
        int atomicN = getAtomicNumber(atomSym);
        Ne = setNumberOfElements(rN,2, order); 
        //Creates the asfem object
        ASFPM asfem{Ne,order,meshType,atomicN,charge,angular,rN};
    
        asfem.setData(confType);
        // **** STARTS THE PROGRAM:  
        // 1. call the load() function:
        // 2. call the size method to get the size of the FEM method
        asfem.load(); 
        int bcDomain = asfem.getBCSize();
        int allOcc = atomicN/2;
        // The SCF Method needs the core matrix hij, build it by calling a first diagonalization:
        // Don't forget to allocate the matrix
    
        //Matrices allocation 
        double *hij = new double[bcDomain*bcDomain];
        double *densMat = new double[bcDomain*bcDomain];
        double *fij = new double[bcDomain*bcDomain];
        double *rho = new double[bcDomain];
        double *rhox = new double[bcDomain];
        double *hpot = new double[bcDomain];
        double *xpot = new double[bcDomain];
        double numElectrons = 4.0;
        double energyHF;
        //Exchange and Hartree potentials: 
        double *vh{nullptr}, *vx{nullptr};

        //Strongly reccomended to fill with zeroes the matrices
        FillZeroMat(fij,bcDomain,bcDomain);
        FillZeroMat(hij,bcDomain,bcDomain);
        FillZeroMat(densMat,bcDomain,bcDomain);



        //Build the core matrix by performing a single diagonalization:
        asfem.singleDiagonalization(hij); // This method gives you the WFN and The orbital Energies
        
        asfem.wfnNormalization(); 
        //
        double *orbitals = asfem.getOrbitals(allOcc);
        //asfpm.printWfn(1);
       

        //Performs the first approach to the Hartree-Fock energy with no electron-electron terms


        asfem.closedShell::computeDensityMatrix(densMat,orbitals);// Builds the density matrix
        energyHF = asfem.closedShell::solveHartreeFock(fij,hij,densMat);
        printf("First Hartree-Fock energy = %lf\n",energyHF);
        //Generates a firtst approach to the vh: Hartree Potential and the vx: exchange potential

        vh = asfem.computeHartreePotential(); 
        vx = asfem.computeExchangePotential();

        for (int i = 0; i < asfem.getDomainSize(); i++)
        {
            printf("vh = % lf\t vx = % lf\n",vh[i], vx[i]);
        }
        asfem.sfc(fij, hij, vh, vx, energyHF);


        /* asfem.getDensity(rho);
        asfem.getExchangeDens(rhox,0,1);

        asfem.solvePoissonEquation(hpot,rho,numElectrons);
        asfem.solvePoissonEquation(xpot,rhox,0.f);

        

        for(int i=0; i<asfem.getBCSize(); i++){
            double exdens = orbitals[i + 0*bcDomain]*orbitals[i + 1*bcDomain];
            double dens1 = orbitals[i + 0*bcDomain]*orbitals[i + 0*bcDomain];
            double dens2 = orbitals[i + 1*bcDomain]*orbitals[i + 1*bcDomain];
            printf("R*rho = % lf    hpot = % lf      R*rhox = % lf     xpot = % lf\n",rho[i],hpot[i],rhox[i],xpot[i]);
        } 
        double mult = 0.0*1000000.0;
        printf("mult = %lf\n",mult); */

        // The SCF Method integrate the potentials above behind the scenes:
        
        //asfpm.sfc(fij, hij, vh, vx, energy);

        //std::cout<<"Final Hartree-Fock energy: "<<energy<<std::endl;
        //std::cout<<std::fixed<<std::setprecision(8)<<"Final Hartree-Fock energy: "<<energy<<std::endl;
        printUsrData(femModel,Ne, order, atomicModel,confType,rC,wallValue,lambda,charge,atomSym,multi,angular,r0,rN,meshType,integrals);
        asfem.printWfn(1);

        delete []densMat;
        delete [] orbitals;
        delete [] fij;
        delete [] hij;
        delete [] rho;
        delete [] hpot;
        delete [] rhox;
        delete [] xpot;
        
    }
    return 0;
}