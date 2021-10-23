#include "closedShell.hpp"

    ClosedShell::ElectronicStructure::ElectronicStructure(){
        std::cout<<"Working with Closed-Shell computations";
    }
    ClosedShell::ElectronicStructure::ElectronicStructure(int _Ne, int _order,Grid<double> _grid, int _atomicN, int _numElec,int _virtualOrbs, double _rMax)
    : Atomic{_Ne, _order, _grid, _atomicN, _numElec, _numElec/2,_virtualOrbs, _rMax}{
       
    }
    ClosedShell::ElectronicStructure::~ElectronicStructure(){
        printf("It's donde!\n");
    }
    //---- Methods ------------
    void ClosedShell::ElectronicStructure::getDensityMatrix(double *densMat){ //May be implemented in the Atomic Structure Static Class
        int k=0;
        for(int i=0; i<bcSize; i++)
        {
            for(int j=0; j<bcSize; j++)
            {   densMat[k] = 0.0;
                for(int orb=0; orb<numOrb; orb++)
                {
                    densMat[k] += 2.0*(orbVec[i + orb*bcSize]*orbVec[j + orb*bcSize]);
                }
                //densMat[k] = 2.0*densMat[k];
                k++;
            }
        }

    }
    void ClosedShell::ElectronicStructure::getTotalDensity(double *rhoinput){
        for(int i=0; i<bcSize; i++){
            rhoinput[i] = 0.0;
            for(int j=0; j<numOrb; j++){
                rhoinput[i] += (orbVec[i + j*bcSize]*orbVec[i + j*bcSize]);  
            }
            rhoinput[i] = 2.0*rhoinput[i];
        }
    }
    void ClosedShell::ElectronicStructure::getPairDensity(double *rhox,int a_orb, int i_orb){
        for(int i=0; i<bcSize; i++){
            rhox[i] = 0.0;
            rhox[i] = rhox[i] + orbVec[a_orb*bcSize + i]*orbVec[i_orb*bcSize + i];
            //printf("pairdens[%d] = %lf\n",i,rhox[i]);
        }
      
    }
    double *ClosedShell::ElectronicStructure::computeHatreePotential(double *uhpot){
        double *vhpot = new double[total_nodes];
    
        for(int i=0; i<bcSize; i++){
            vhpot[i+1]= uhpot[i]/femgrid[i+1];
        }
        vhpot[0]=doInterpolation(&femgrid[0],vhpot,2);
        vhpot[total_nodes-1]=totQ/femgrid[total_nodes-1];
  
        return vhpot;
    }
    double *ClosedShell::ElectronicStructure::computeAuxiliarExchangePotential(int inputOrbital, double *sij){
        double *wx = new double[numOrb*bcSize];
        double *sourceVec = new double[bcSize];
        double *aux = new double[bcSize];
        double *rhoX = new double[bcSize];
        double *rightVector;
        double bundaryConValue;
        for(int a=0; a<numOrb; a++){
            printf("Solving Poisson equation for orbital %d\n",a+1);
            if(a==inputOrbital){
               bundaryConValue = 1.0;
                
            }
            else{
                 bundaryConValue = 0.0;
            }
            printf("Boundary value = %lf\n",bundaryConValue);
            getPairDensity(rhoX,inputOrbital,a);
            rightVector = divideOverGridPoints(rhoX);
            getSourceVector(sourceVec,sij,rightVector,bundaryConValue);
            solvePoissonEquation(aux,sourceVec);
            for(int i=0; i<bcSize; i++){
                wx[i + a*bcSize] = aux[i];
                aux[i]=0.0;
            }
        }
        delete [] sourceVec;
        delete [] aux;
        delete [] rhoX;
        delete [] rightVector;
        return wx;
    }
    void ClosedShell::ElectronicStructure::eigenSystemSCF(double *hcore,double *sij, double *matCoeffs,double *eigenVal, bool flag,int inputTol){
        //Computes a full SCF by solving the EigenSystem Fc=SC; 
        // The matrix F does not contain the Exchange Term
        if(flag==true){
            int tol = 100;
            if(inputTol!=0){
                tol = inputTol;
            }
            //Performs a full SCF
            std::cout<<"Performing a full eigenvalue problem\n";
            double *densmat = new double[bcSize*bcSize];
            double *uhpot = new double[bcSize];
            double *sourceVec = new double[bcSize];
            double *vh_mat = new double[bcSize*bcSize];
            double *fock_mat = new double[bcSize*bcSize];
            double *density = new double[bcSize];
            double *vh,*rightVector;
            diag(bcSize,hcore,sij,eigenVal,matCoeffs);
            wfnNormalization(matCoeffs);
            getOrbitals(matCoeffs); //Obtains the occupied Orbitals
            getTotalDensity(density);
            rightVector = divideOverGridPoints(density);
            getSourceVector(sourceVec,sij,rightVector,totQ);
            solvePoissonEquation(uhpot,sourceVec);
            vh = computeHatreePotential(uhpot);
            double energy0  = 0.5*eigenVal[0];
            int iter = 0;
            double deltaEnergy = 100000000.0;
            double new_HF = energy0;
            double old_HF,orb_energy;
            printf("energy0 = %lf\n",energy0);
            while(deltaEnergy>0.00000001 || iter<tol){
                old_HF = new_HF;
                femNumericalIntegration(vh_mat,vh);
                SumMatrices(hcore,vh_mat,fock_mat,bcSize*bcSize);
                diag(bcSize,fock_mat,sij,eigenVal,matCoeffs);
                orb_energy = 0.5*eigenVal[0];
                wfnNormalization(matCoeffs);
                getOrbitals(matCoeffs);
                getDensityMatrix(densmat);
                new_HF = energyHF(hcore,fock_mat,densmat);
                printf("Hartree-Fock energy = %.10lf    orb = %.10lf   iteration: %d\n",new_HF,orb_energy,iter+1);
                getTotalDensity(density);
                rightVector = divideOverGridPoints(density);
                getSourceVector(sourceVec,sij,rightVector,totQ);
                solvePoissonEquation(uhpot,sourceVec);
                vh  = computeHatreePotential(uhpot);
                deltaEnergy = fabs(new_HF-old_HF);
                iter++;
            }
            samePhase(0,1);
            printWfn();
            /*femNumericalIntegration(vh_mat,vh); */
            //vh = computeHatreePotential()


            delete [] density;
            delete [] densmat;
            delete [] vh;
            delete [] rightVector;
            delete [] sourceVec;
            delete [] uhpot;
            delete [] vh_mat;
        }
        else{
            double energy0;
            diag(bcSize,hcore,sij,eigenVal,matCoeffs);
            printf("energy0 = %.16lf\n",0.5*eigenVal[0]);
            //asfem_tools::diag(matSize,hcore,sij,eigenVal,wfn);
        }
        
    }
    void ClosedShell::ElectronicStructure::iterativeSCF(double *hcore, double *sij, double *matCoeffs){
        //Perform a few iterations solving an eigenvalue problem:
        double *uhpot = new double[bcSize];
        double *sourceVec = new double[bcSize];
        double *vh_mat = new double[bcSize*bcSize];
        double *fock_mat = new double[bcSize*bcSize];
        double *density = new double[bcSize];
        double *eigenVal = new double[bcSize];
        double *vh,*rightVector;
        double *wx,*vx,*exchangeVec;
        diag(bcSize,hcore,sij,eigenVal,matCoeffs);
        int tol = 0; 
        while(tol<5){
            printf("Diagonalizatoin step %d\n",tol);
            wfnNormalization(matCoeffs);
            getOrbitals(matCoeffs); //Obtains the occupied Orbitals
            getTotalDensity(density);
            rightVector = divideOverGridPoints(density);
            getSourceVector(sourceVec,sij,rightVector,totQ);
            solvePoissonEquation(uhpot,sourceVec);
            vh = computeHatreePotential(uhpot);
            femNumericalIntegration(vh_mat,vh);
            SumMatrices(hcore,vh_mat,fock_mat,bcSize*bcSize);
            diag(bcSize,fock_mat,sij,eigenVal,matCoeffs);
            tol++;
        }
        printf("Eigenvalue by diagonalization  =  %lf\n",eigenVal[0]);
        wx = computeAuxiliarExchangePotential(0,sij);
        vx = computeExchangePotential(wx,0);
        exchangeVec = integrateExchangePotential(vx);
        rayleighQuotient(fock_mat,sij,exchangeVec);
        /* for(int i=0; i<bcSize; i++){
            printf("exVec[%d] = %lf\n",i, exchangeVec[i]);
            //printf("wx12 = %lf     wx22  = %lf\n",wx[i + 0*bcSize], wx[i + 1*bcSize]);
        } */
        //getOrbitals(matCoeffs);  //Gets the newest orbitals
        /* vx = computeExchangePotential(wx,1);
        for(int i=0; i<globalSize; i++){
            printf("r = % lf    VH = % lf     Vx21[%d] =  % lf    Vx22[%d] = % lf\n",femgrid[i],vh[i],i, vx[i + 0*globalSize],i,vx[i + 1*globalSize]);
        } */
        //Starting the iterative procedure: 
        samePhase(0,1);
        printWfn();   
        //rayleighQuotient(hcore,sij,vh);
        delete [] density;
        delete [] wx;
        delete [] vh;
        delete [] rightVector;
        delete [] sourceVec;
        delete [] uhpot;
        delete [] vh_mat;

    }
