#include "atomicStructure.hpp"


Atomic::Atomic(){
    //std::cout<<"Atomic Structure Computations\n";
}
Atomic::Atomic(int _Ne, int _order,Grid<double> _grid, int _atomicN, int _numElec, int _orbs,int _virtualOrbs,double _rMax)
: atomicN{_atomicN},numOrb{_orbs}, r0{0}, rN{_rMax},numElec{_numElec},virtualOrbs{_virtualOrbs}, FixedPoints{_Ne,_order,_grid}{
    std::cout<<"Atomic Structure Computations\n";
    total_nodes = Ne*order+1;
    fem_nodes = total_nodes - 2;
    totQ = static_cast<double>(_numElec); 
    orbVec = new double[(numOrb + virtualOrbs)*fem_nodes];
    
}
Atomic::~Atomic(){/*Destructor*/
    printf("Have a nice Day\n");
    delete [] orbVec;
}
void Atomic::getOrbitals(double *matCoeff){
    int k=0;
    for(int i=0; i<(numOrb+virtualOrbs); i++){
        for(int j=0; j<bcSize; j++){
            orbVec[k] = matCoeff[i*bcSize + j];
            k++;
        }
    }
}
int Atomic::orbitalPhase(int orbOfinterest){
    int i, answer,phase; 
    double prev, val, der;
    prev = orbVec[orbOfinterest*bcSize];
    i=1;
    answer = 1;
    while(i<bcSize && answer!=0){
        val = orbVec[i + orbOfinterest*bcSize];
        der = val-prev;
        if(fabs(der)>1E-5){
            answer=0;
            //printf("The orbital changes the phase in %d\n",i);
        }
        prev = val;
        i++;
    }
    if(der>0){
        return phase = 1;
    }
    else{
        return phase = -1;
    }
    return 0;
}
void Atomic::printOrbital(std::string name, int target){
    std::cout<<"Orbital "<<name<<" data: \n";
    for(int i=0; i<bcSize; i++){
        //std::cout<<name<<"["<<i<<"]"<<orbVec[i + target*bcSize]<<std::endl;
        printf("orb[%d] = %lf\n",i,orbVec[i + target*bcSize]);
    }
}
void Atomic::printWfn(){
    //int realOrb = orbital-1;
    double wfnOrb;
    int phase;
    int phase2;
    std::ofstream wfnData;
    wfnData.open("wfn.dat",std::ios::out);
    if(!wfnData){
        std::cout<<"FILE NOT CREATED\n";
    } 
    phase = orbitalPhase(0);
    phase2 = orbitalPhase(1);
    std::cout<<"The orbital 1s Phase is: "<<phase<<std::endl;
    std::cout<<"The orbital 2s Phase is: "<<phase2<<std::endl;
    wfnData<<"#r \t 1s \t 2s "<<std::endl;
    wfnData<<std::fixed<<std::setprecision(6)<<femgrid[0]<<"\t"<<0.0<<"\t"<<0.0<<std::endl;
    for(int i=0; i<bcSize; i++)
    {
        wfnData<<std::fixed<<std::setprecision(6)<<femgrid[i+1]<<"\t"<<orbVec[i + bcSize*0]<<"\t"<<orbVec[i + bcSize*1]<<std::endl;

    }
    //wfnData<<std::fixed<<std::setprecision(6)<<femgrid[0]<<"\t"<<0.0<<"\t"<<0.0<<std::endl;
    wfnData.close();
}
double Atomic::energyHF(double *hij, double *fij,double *densMat){
    
    double energy = 0.0;
    for(int i=0; i<bcSize*bcSize; i++)
    {
        energy = energy + 0.5*(densMat[i]*hij[i]);
        energy = energy + 0.5*(fij[i]*densMat[i]);
    }
    return 0.5*energy;

}
double *Atomic::divideOverGridPoints(double *input){
    double *output = new double[bcSize];
    for(int i=0; i<bcSize; i++){
        output[i] = input[i]/femgrid[i+1];
    }
    return output;
}
// *** Methods ***
void Atomic::wfnNormalization(double *wfn){
    //Normalizes the full wavefunction
    double *cf = new double[total_nodes];
    double *feMatS;
    double normConst;
    double x[poly];
    for(int orbital=0; orbital<fem_nodes; orbital++){
        cf[0] = 0.0; //Boundary condition 1
        cf[total_nodes-1] = 0.0; //Boundary Condition 2
        
        for(int i=1; i<fem_nodes+1; i++){
            cf[i] = wfn[(i-1) + fem_nodes*orbital];
            //printf("cf = %lf\n",cf[i]);
        }
        normConst=0.0;
        for(int ei=0; ei<Ne; ei++){
            for(int j=0; j<poly; j++){
                int indx = poly*ei + j;
                int i = linkMat[indx];
                x[j] = femgrid[i];
            }
            feMatS = getFixedPointsOverlapMatrices(x,order);
            normConst  += integrateElement(ei,order,feMatS,&linkMat[0],1.0, cf);
        }
        normConst = 1.0/sqrt(normConst);
        //printf("consN = %lf\n",normConst);
        for(int i=0; i<fem_nodes; i++){
            wfn[i + fem_nodes*orbital] *= normConst;
            //printf("Wfn = %lf\n",wfn[i+fem_nodes*orbital]);
        } 
    }
    delete [] feMatS;
    delete [] cf;
}
void Atomic::rayleighQuotient(double *f_mat, double *smat, double *exchangevec){
    int totOrbitals = (numOrb+virtualOrbs); //Number of total orbitals selected by the user
    double *trialVec = new double[(numOrb+virtualOrbs)*bcSize];
    double *matprodNum= new double[(numOrb+virtualOrbs)*bcSize];
    double *matprodDenom= new double[(numOrb+virtualOrbs)*bcSize];
    double *orbitalEnergy = new double[(numOrb+virtualOrbs)];
    for(int i=0; i<(numOrb + virtualOrbs); i++){
        for(int j=0; j<bcSize; j++){
            trialVec[j + i*bcSize] = orbVec[j + i*bcSize]; //Takes the trial vector;
        }
    }
    matMult(trialVec,totOrbitals,bcSize,f_mat,bcSize,bcSize,matprodNum);
    matMult(trialVec,totOrbitals,bcSize,smat,bcSize,bcSize,matprodDenom);
    double numerator,denominator;
    for(int i=0; i<(numOrb + virtualOrbs); i++){
        numerator = 0.0;
        denominator = 0.0;
        for(int j=0; j<bcSize; j++){
            
            //auxnum[j + i*bcSize] = auxnum[j + i*bcSize];
            //printf("vx = %lf\n",exchangevec[j]);
            numerator = numerator + trialVec[j + i*bcSize]*matprodNum[j + i*bcSize];
            denominator = denominator + trialVec[j + i*bcSize]*matprodDenom[j + i*bcSize];
            
        }
        orbitalEnergy[i] = numerator/denominator;
        printf("Orbital %d energy iterative SCF = %lf\n",i,0.5*orbitalEnergy[i]);
    } 
    delete [] trialVec;
    delete [] orbitalEnergy;
    delete [] matprodNum;
     delete [] matprodDenom;
}
void Atomic::samePhase(int orb1, int orb2){
    int phase1 = orbitalPhase(orb1);
    int phase2 = orbitalPhase(orb2);

    for(int i=0; i<bcSize; i++){
        orbVec[i + orb1*bcSize] = phase1*orbVec[i + orb1*bcSize];
        orbVec[i + orb2*bcSize] = phase2*orbVec[i + orb2*bcSize];
    }
}
double * Atomic::computeExchangePotential(double *wx, int inputOrb){
    double *vx = new double[numOrb*globalSize];
    double *aux = new double[globalSize];
    aux[0]=0.0; aux[globalSize-1] = 0.0;
    for(int a=0; a<numOrb; a++){
        for(int i=0; i<bcSize; i++){
            aux[i+1] = 0.0;
            aux[i+1] = wx[i + a*globalSize]/femgrid[i+1];
        }
        aux[0]=doInterpolation(&femgrid[0],aux,2);
        if(a!=inputOrb){
            aux[total_nodes-1] = 0.0;
        }
        else{
            aux[total_nodes-1] = 1.0/femgrid[total_nodes-1];
        }
        for(int i = 0; i<globalSize; i++){
            vx[i + a*globalSize] = aux[i];
        }
    }
    
    delete [] aux;
    return vx;
}
double * Atomic::integrateExchangePotential(double *wx){
    double *vx = new double[bcSize];
    double *aux_vec = new double[globalSize];
    double *matEx = new double[bcSize*bcSize];
    double *vec = new double [bcSize];
    double *vecprod = new double [bcSize];
    FillZeroMat(matEx,bcSize,bcSize);
    for(int a=0; a<numOrb; a++){
        //Performs an integration for each orbital
        
        for(int i=0; i<globalSize; i++){
            aux_vec[i]  = wx[i + a*globalSize];
        }
        femNumericalIntegration(matEx,aux_vec);
        for(int i=0; i<bcSize; i++){
            vec[i] = orbVec[i + a*bcSize];
            vecprod[i]=0.0;
        }
        MatrixProduct(matEx,vec,vecprod,bcSize,bcSize,1);
        //printf("Multiplication of orbital %d \n", a);
        for (int i = 0; i < bcSize; i++)
        {
            //printf("prod = %lf\n", vecprod[i]);
            vx[i] = vx[i] + vecprod[i];
        }
        
    }
    delete [] aux_vec;
    delete [] vecprod;
    delete [] vec;
    delete [] matEx;
    return vx;
}