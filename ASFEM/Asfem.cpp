#include<iostream>
#include "Asfem.hpp"
#include<string>
// External Potential Integration options
// "Numeric" it needs the External Potential Evaluated in the FEM GRID
// Numeric integration needs an extrapolated value at V(r=0)
// "Analitic" Computes the numeric integrals avoiding the V(R=0)
// ** CONSTRUCTORS *****
ASFEM::ASFEM(){/*Default Constructor*/}

ASFEM::ASFEM(std::string femModel, int Ne, int order,std::string _atomicModel, double _lambda,std::string _confType,double _Rc,double _wallVal,int _atomicN, int _charge,int _angular, std::string gridType,double rInfty,std::string _integrals)
:FEM{Ne,order,femModel, gridType},r0{0.0},rN{rInfty},atomicN{_atomicN},atomicModel{_atomicModel},integrationScheme{_integrals},wallValue{_wallVal},
cutRad{_Rc},lambda{_lambda},charge{_charge},angular{_angular}{
    std::cout<<"ASFEM constructor 1.\n";
    total_nodes = Ne*order+1;
    fem_nodes = total_nodes-2;
    numElectrons = atomicN-charge;
    totQ = numElectrons;
    occOrb = atomicN/2;
    wfn = new double[bcSize*bcSize];
    eigenVal = new double[bcSize];
    rho = new double[bcSize];
    vr = new double[globalSize];

}

ASFEM::ASFEM(int Ne, int order, int poissonNe, 
std::string femModel,double _r0, 
double _rN,std::string gridType,int atom, std::string intAtomicModel,std::string potInteg)
:FEM{Ne, order, poissonNe,femModel,gridType},r0{_r0},rN{_rN},atomicN{atom},atomicModel{intAtomicModel},integrationScheme{potInteg}{
    total_nodes = Ne*order +1;
    fem_nodes = total_nodes - 2;
    occOrb = atomicN/2;
    numElectrons = atomicN-charge;
    totQ = numElectrons;
    vr = new double[globalSize];
    wfn = new double[bcSize*bcSize];
    eigenVal = new double[bcSize];
    rho = new double[bcSize];
}
ASFEM::~ASFEM(){

    delete [] vr;
    delete [] wfn;
    delete [] eigenVal;
    delete [] rho;
}
// *** END CONSTRUCTORS AND DESTRUCTORS
//******** PRIVATE METHODS ******************
void ASFEM::wfnNormalization(std::string name){
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
        /* printf("consN = %lf\n",normConst); */
        for(int i=0; i<fem_nodes; i++){
            wfn[i+fem_nodes*orbital] *= normConst;
            //printf("Wfn = %lf\n",wfn[i+fem_nodes*orbital]);
        } 
    }
    delete [] cf;
    delete [] feMatS;
}
void ASFEM::wfnNormalization(){
    double *cf = new double[total_nodes];
    double *feMatS;
    feMatS = ElementalOverlapMatrix(order);
    double normConst;
    for(int orbital=0; orbital<fem_nodes; orbital++){
        cf[0] = 0.0; //Boundary condition 1
        cf[total_nodes-1] = 0.0; //Boundary Condition 2
        
        for(int i=1; i<fem_nodes+1; i++){
            cf[i] = wfn[(i-1) + fem_nodes*orbital];
            //printf("cf = %lf\n",cf[i]);
        }
        normConst=0.0;
        for(int ei=0; ei<Ne; ei++){
            double coeff = 0.5*femgrid.getElementSize(ei);
            normConst  += integrateElement(ei,order,feMatS,&linkMat[0],coeff, cf);
        }
        normConst = 1.0/sqrt(normConst);
        for(int i=0; i<fem_nodes; i++){
            wfn[i+fem_nodes*orbital] *= normConst;
            //printf("Wfn = %lf\n",wfn[i+fem_nodes*orbital]);
        } 
    }
    delete [] cf;
    delete [] feMatS;
}
// ***** COMPUTE DENSITY AND DENSITY MATRIX**************************************
void ASFEM::getDensity(){
    if(numElectrons!=1){
        for(int i=0; i<bcSize; i++){
            for(int orb=0; orb<occOrb; orb++){
                rho[i] += wfn[i + orb*bcSize]*wfn[i + orb*bcSize];
                rho[i] = 2.0*rho[i];
            }
        }
    }
    else{//This is the density for the Hidrogen like Atoms;
        int orb=0;
        for(int i=0; i<bcSize; i++){
            rho[i] = wfn[i + orb*bcSize]*wfn[i+orb*bcSize];
        }
    }
}
void ASFEM::getDensityMatrix(double *densMat){ //May be implemented in the Atomic Structure Static Class
    int k=0;
    std::cout<<"Orbital od Interest: "<<occOrb<<std::endl;
    for(int i=0; i<bcSize; i++){
        for(int j=0; j<bcSize; j++){
            for(int orb=0; orb<occOrb; orb++){
                densMat[k] = 2.0*(wfn[i + orb*bcSize]*wfn[j + orb*bcSize]);
            }
            k++;
        }
    }

}
//**************************************************************************************
//***** HARTREE_POTENTIAL **********************
double * ASFEM::computeHartreePotential(double *hpot){
    double *vh_r = new double[total_nodes];

    for(int i=0; i<bcSize; i++){
        vh_r[i+1]= hpot[i]/femgrid[i+1];
    }
    vh_r[0]=doInterpolation(&femgrid[0],vh_r,4);
    vh_r[total_nodes-1]=totQ/femgrid[total_nodes-1];
    return vh_r;
}
//*******************************************
double ASFEM::energyHF(double *hij, double *fij,double *densMat){
    
    double energy = 0.0;
	for(int i=0; i<bcSize*bcSize; i++)
	{
		energy = energy + 0.5*(hij[i]*densMat[i]);
		energy = energy + 0.5*(fij[i]*densMat[i]);
	}
	return 0.5*energy;

}
void ASFEM::singleDiagonalization(){
    double *hij = new double[bcSize*bcSize];
    int phase;
    SumMatrices(&vij[0],&kij[0],hij,bcSize*bcSize);
    diag(bcSize,hij,&sij[0],eigenVal,wfn);
    //asfem_tools::sumArrays(hij,&vij[0],&kij[0],bcSize,bcSize);
    //asfem_tools::diag(bcSize,hij,&sij[0],eigenVal,wfn);
    printf("orbital value: %.10lf\n",0.5*eigenVal[0]);
    wfnNormalization(femModel);
    asfem_tools::wfn::getWfnPhase(fem_nodes,0,&phase,wfn);
    delete [] hij;
}
void ASFEM::getExternalPotential(){

    evaluateExternalPotential(vr,&femgrid[0], 0, atomicN, globalSize);
    vr[0] = asfem_tools::doInterpolation(&femgrid[0],vr,4);
    //vr[0] = 14.447169; 
    printf("Extrapolation Value v(r=0) = %lf\n",vr[0]);
}
//************* SCF*************************************
void ASFEM::performSCF(){
    double *hij = new double[bcSize*bcSize];
    double *vhij = new double[bcSize*bcSize];
    double *densMat = new double[bcSize*bcSize];
    double *fij = new double[bcSize*bcSize];
    double *R_rho = new double[bcSize];
    double *hpot = new double[bcSize];
    double *R_hpot{nullptr};
    FillZeroMat(fij,bcSize,bcSize);
    FillZeroMat(vhij,bcSize,bcSize);
    FillZeroMat(hij,bcSize,bcSize);
    int phase;
    SumMatrices(&vij[0],&kij[0],hij,bcSize*bcSize);
    diag(bcSize,hij,&sij[0],eigenVal,wfn);
    printf("First orbital value: %.10lf\n",0.5*eigenVal[0]);
    wfnNormalization(femModel);
    asfem_tools::wfn::getWfnPhase(fem_nodes,0,&phase,wfn);
    getDensityMatrix(densMat);
    double energy0 = energyHF(hij,fij,densMat);
    printf("First Hartree-Fock Energy: %.10lf\n",energy0);
    getDensity();
    divideBy(R_rho,rho,&femgrid[0],bcSize);
    solvePoissonEquation(hpot, R_rho,totQ);
    R_hpot = computeHartreePotential(hpot);
    double new_hf= energy0;
    double orb_e,old_hf;
    //fixedElementsNumIntegration(vhij,R_hpot);
    fixedPointsNumIntegration(vhij,R_hpot);
    /* for(int i=0; i<bcSize; i++){
        printf("dens = %lf\n",rho[i]);
    } */

    SumMatrices(hij,vhij,fij,bcSize*bcSize);
    diag(bcSize,fij,&sij[0],eigenVal,wfn);
    printf("Orbital value: %.10lf\n",0.5*eigenVal[0]);

    delete [] hij;
    delete [] vhij;
    delete [] densMat;
    delete [] fij;
    delete [] R_rho;
    delete [] hpot;
    delete [] R_hpot;
    

    
}
//****** SCF**********************************************
//********  PUBLIC METHODS **************************
void ASFEM::printWfn(int orbital){
    int realOrb = orbital-1;
    double wfnOrb;
    int phase;
    std::ofstream wfnData;
    wfnData.open("wfn.dat",std::ios::out);
    if(!wfnData){
        std::cout<<"FILE NOT CREATED\n";
    } 
    asfem_tools::wfn::getWfnPhase(fem_nodes,realOrb,&phase,wfn);
    std::cout<<"The Wave Function Phase is: "<<phase<<std::endl;
    for(int i=0; i<fem_nodes; i++)
    {
        wfnData<<std::fixed<<std::setprecision(6)<<wfn[i + fem_nodes*realOrb]*(phase)<<std::endl;

    }
    wfnData.close();
}
void ASFEM::startProgram(){
    //This public Method starts the ASFEM program
    if(atomicN==1){
        singleDiagonalization();
    }
    //Build the FEM Grid
     if(femModel=="Fixed Elements"){
        buildFemGrid(atomicN,r0,rN);
        getExternalPotential();
        assambleMatricesFixedElements(vr);
        //singleDiagonalization();
         performSCF();
    }
    else{
        buildFemGrid(atomicN,r0,rN);
        assambleMatricesFixedPoints(atomicN);
        performSCF();
    }
}
void ASFEM::printInitialData(){

    std::cout<<"********** ASFEM INITIAL DATA **********\n";
    std::cout<<"Finite Element Model: "<<femModel<<std::endl;
    std::cout<<"Number of Elements: "<<Ne<<std::endl;
    std::cout<<"Total Points: "<<(Ne*order + 1)<<std::endl;
    std::cout<<"Polynomial Order: "<<order<<std::endl;
    std::cout<<"Atomic Number: "<<atomicN<<std::endl;
    std::cout<<"Model: "<<atomicModel<<std::endl;
    if(atomicModel=="Confined-Atom"){
        std::cout<<"Confinement: "<<confType<<std::endl;
        if(confType=="Soft-Walls"){
            std::cout<<"Cut Radii: "<<cutRad<<std::endl;
            std::cout<<"Wall-Value: "<<wallValue<<std::endl;
        }
        else{
            std::cout<<"Cut Radii: "<<cutRad<<std::endl;
        }
    }
    else if(atomicModel=="Plasma"){
        std::cout<<"Lambda Value: "<<lambda<<std::endl;
    }
    std::cout<<"Practical Infinity: "<<rN<<std::endl;
    std::cout<<"Grid Form: "<<gridType<<std::endl;
    std::cout<<"External Potential Integral: "<<integrationScheme<<std::endl;
    std::cout<<"***************************************\n";
}