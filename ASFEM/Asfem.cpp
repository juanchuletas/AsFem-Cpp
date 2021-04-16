#include<iostream>
#include "Asfem.hpp"
#include<string>
// External Potential Integration options
// "Numeric" it needs the External Potential Evaluated in the FEM GRID
// Numeric integration needs an extrapolated value at V(r=0)
// "Analitic" Computes the numeric integrals avoiding the V(R=0)
ASFEM::ASFEM(){/*Default Constructor*/}

ASFEM::ASFEM(int Ne, int order, int poissonNe, 
std::string femModel,double _r0, 
double _rN,std::string gridType,int atom, std::string intAtomicModel,std::string potInteg)
:FEM{Ne, order,poissonNe, femModel,gridType},r0{_r0},rN{_rN},atomicN{atom},atomicModel{intAtomicModel},integrationScheme{potInteg}{
    total_nodes = Ne*order +1;
    fem_nodes = total_nodes - 2;
    occOrb = atomicN/2;
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
            feMatS = getOverlapElementalMatrices(x,order);
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
            printf("coef = %lf\n",coeff);
            double phi[poly];
            for(int j=0; j<poly; j++){
                phi[j] = cf[linkMat[poly*ei+j]];
                /* printf("phi = %lf\n",phi[j]); */
            }
            double value = 0.0;
            for(int mu=0; mu<poly; mu++){
                for(int nu=0; nu<poly; nu++){
                    value = value + phi[mu]*phi[nu]*(coeff*feMatS[poly*mu + nu]);
                    //printf("Value = %lf\n",normConst);
                }
            }
            normConst += value;
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
void ASFEM::getDensity(){
    if(atomicN!=1){
        for(int i=0; i<bcSize; i++){
            for(int orb=0; orb<occOrb; orb++){
                rho[i] += wfn[i + orb*bcSize]*wfn[i + orb*bcSize];
                rho[i] += 2.0*wfn[i];
            }
        }
    }
    else{//This is the density for the Hidrogen Atom;
        int orb=0;
        for(int i=0; i<bcSize; i++){
            rho[i] = wfn[i + orb*bcSize]*wfn[i+orb*bcSize];
        }
    }
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
    for(int i=0; i<bcSize; i++){
        printf("WFN NORMALIZED = %lf\n",wfn[i + bcSize*0]*phase);
    } 
    delete [] hij;

}
void ASFEM::getExternalPotential(){

    evaluateExternalPotential(vr,&femgrid[0], 0, atomicN, globalSize);
    vr[0] = asfem_tools::doInterpolation(&femgrid[0],vr,4);
    //vr[0] = 14.447169; 
    printf("Extrapolation Value v(r=0) = %lf\n",vr[0]);
}
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

    //Build the FEM Grid
    if(femModel=="Fixed Elements"){
        buildFemGrid(atomicN,r0,rN);
        getExternalPotential();
        singleDiagonalization();
    }
    else{
        buildFemGrid(atomicN,r0,rN);
        assambleFemMatrices(atomicN);
        singleDiagonalization();
    }
}
void ASFEM::printInitialData(){

    std::cout<<"********** ASFEM INITIAL DATA **********\n";
    std::cout<<"Finite Element Model: "<<femModel<<std::endl;
    std::cout<<"Number of Elements: "<<Ne<<std::endl;
    std::cout<<"Total Points: "<<(Ne*order + 1)<<std::endl;
    std::cout<<"Polynomial Order: "<<order<<std::endl;
    std::cout<<"Atom: "<<atomName<<std::endl;
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