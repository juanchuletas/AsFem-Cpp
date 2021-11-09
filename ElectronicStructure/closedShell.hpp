#if !defined(_ElectronicStructure_H_)
#define _ElectronicStructure_H_
#include <iostream>
#include<iostream>
#include "../Atomic_Computations/atomicStructure.hpp"
namespace ClosedShell{

    class ElectronicStructure : public Atomic{
            /* int occOrbs;
            int numElectrons;
            int atomicN;
            int gridSize;
            int bcDomSize;
            int matSize; */
            /* double *orbital;
            double *orbitalEnergy; */
        public:
            //** Constructors & Destructors
            ElectronicStructure();
            ElectronicStructure(int _Ne, int _order,Grid<double> _grid, int atomicN, int _numElec,int _virtualOrbs ,double _rMax);
            ~ElectronicStructure();
            // ****  Methods ***
            void eigenSystemSCF(double *hcore,double *sij, double *matCoeffs,double *eigenVal, bool flag,int inputTol);
            void iterativeSCF(double *hcore, double *sij, double *matCoefss);
            void solveHatreeFockequation(double *orbital);
            double * computeHatreePotential(double *iput);
            void getPairDensity(double *xdens,int a_orb, int i_orb);
            void getTotalDensity(double *rhoinput);
            void getDensityMatrix(double *densMat);
            double *computeAuxiliarExchangePotential(int inputOrbital,double *sij);
    };
}


#endif // _ElectronicStructure_H_
