#if !defined(_ELECTRONIC_STRUCTURE_H_)
#define _ELECTRONIC_STRUCTURE_H_
#include <iostream>
#include<string>

namespace electronicStructure{

    class closedShell {
        
        
        int nucleii;
        int charge;
        int orbitals;
        int numElectrons;
        int gridnodes;
        int bcnodes;
        public: 
        closedShell();
        closedShell(int _nuc, int _charge);
        ~closedShell();
        void singleDiagonalization(double *hij);
        void setNumOfGridPoints(int _gridPoints);
        void computeDensityMatrix(double *densMat,double *wfn);
        double solveHartreeFock(double *fij, double *hij, double *densMat);


    };



}
#endif // _ELECTRONIC_STRUCTURE_H_
