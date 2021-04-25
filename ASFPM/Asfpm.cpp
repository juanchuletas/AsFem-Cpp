#include "Asfpm.hpp"


ASFPM::ASFPM(){

}
ASFPM::ASFPM(std::string femModel, int Ne, int order,std::string _atomicModel, double _lambda,std::string _confType,double _Rc,double _wallVal,int _atomicN, int _charge,int _angular, std::string gridType,double rInfty,std::string _integrals)
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
ASFPM::~ASFPM(){

    delete [] wfn;
    delete [] eigenVal;
    delete [] rho;
    delete [] vr;
}
// ********* PUPUPUPUPUPUPUPUPUPUBLIC METHODSSSS ************
void ASFPM::wfnNormalization(){
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
void ASFPM::getDensity(){
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
void ASFPM::getDensityMatrix(double *densMat){ //May be implemented in the Atomic Structure Static Class
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
double * ASFPM::computeHartreePotential(double *hpot){
    double *vh_r = new double[total_nodes];

    for(int i=0; i<bcSize; i++){
        vh_r[i+1]= hpot[i]/femgrid[i+1];
    }
    vh_r[0]=doInterpolation(&femgrid[0],vh_r,4);
    vh_r[total_nodes-1]=totQ/femgrid[total_nodes-1];
    return vh_r;
}