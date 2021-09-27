#if !defined(_CLOSED_SHELL)
#define _CLOSED_SHELL
#include<iostream>
 class ClosedShell{
        int OccupiedOrb;
        int fock;
        int bcDomSize;
        int fullDomSize;
        //This class or namespace solves the Hartree-Fock equations for closed shell systems
        public:
            ClosedShell(int a);
            ~ClosedShell();
            void solveHartreFockEquations();

};

#endif // _CLOSED_SHELL
