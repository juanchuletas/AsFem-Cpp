#include "Asfpm.hpp"


ASFPM::ASFPM(){
}
ASFPM::ASFPM(int _Ne, int _order,std::string _grid,int _atomicN,int _charge, int _angular,double _rN)// Pretty good constructor
: FEMFP{_Ne,_order,"Fixed-Points",_grid},electronicStructure::closedShell{_atomicN,_charge},r0{0.0},
    atomicN{_atomicN},charge{_charge},angular{_angular},rN{_rN}{
    std::cout<<"ASFEM constructor 1.\n";
    
    total_nodes = globalSize;
    fem_nodes = bcSize;
    wfn = new double[fem_nodes*fem_nodes];
    eigenVal = new double[fem_nodes];
    rho = new double[fem_nodes];
}
ASFPM::~ASFPM(){

    delete [] wfn;
    delete [] eigenVal;
    delete [] rho;
}
// SOME SETTERS AND GETTERS
void ASFPM::setData(std::string _confine,double rc, double wall){//Soft-Walls and hard Walls setter
    confType = _confine;
    cutRad = rc;
    wallValue = wall;
    occOrb = atomicN/2;
    numElectrons = atomicN-charge;
    totQ = numElectrons;
    closedShell::setNumOfGridPoints(total_nodes);
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
        //printf("consN = %lf\n",normConst);
        for(int i=0; i<fem_nodes; i++){
            wfn[i+fem_nodes*orbital] *= normConst;
            //printf("Wfn = %lf\n",wfn[i+fem_nodes*orbital]);
        } 
    }
    delete [] cf;
    delete [] feMatS;
}
double *ASFPM::getOrbitals(int numofOrb){
    
    double *orbitals = new double[numofOrb*bcSize];
    int k=0;
    for(int i=0; i<numofOrb; i++){
        for(int j=0; j<bcSize; j++){
            orbitals[k] = wfn[i*bcSize + j];
            k++;
        }
    }
    return orbitals;
}
double ASFPM::computeTotalCharge(int rcIndex){
    int nele = (rcIndex-1)/order;
    nele = nele + 1;
    double *totrho =  new double[total_nodes];
    totrho[0] = 0.0;
    totrho[total_nodes-1] = 0.0;
    double x[poly],phi[poly];
    double *feMatS;
    for(int i=1; i<fem_nodes+1; i++){
        totrho[i] = wfn[i-1];
        //printf("cf = %lf\n",cf[i]);
    }
    double finalq = 0.0;
    for(int ei=0; ei<nele; ei++){
            //printf("e = %d\n",ei);
            for(int j=0; j<poly; j++){
                int indx = poly*ei + j;
                int i = linkMat[indx];
                x[j] = femgrid[i];
                phi[j] = totrho[i];
            }
            feMatS = getFixedPointsOverlapMatrices(x,order);
            double qtot = 0.0;
            for(int mu=0; mu<poly; mu++){
                for(int nu=0; nu<poly; nu++){
                 qtot = qtot + 2.0*phi[mu]*phi[nu]*(feMatS[poly*mu + nu]);
                 //printf("mats[%d] = %lf\n",poly*mu + nu,feMatS[poly*mu + nu]);
                }
            }
            finalq = finalq + qtot; 
            
    }
    delete [] totrho;
    delete [] feMatS;
    return finalq;
} 
void ASFPM::getDensity(const double occup){

    for(int i=0; i<bcSize; i++){
        rho[i] = 0.0;
        for(int orb=0; orb<occOrb; orb++){
            rho[i] += wfn[i + orb*bcSize]*wfn[i + orb*bcSize];
          
        }
        rho[i] = occup*rho[i];
    }
}
void ASFPM::getDensity(double *rhor){
    for(int i=0; i<bcSize; i++){
        rhor[i] = 0.0;
        for(int orb=0; orb<occOrb; orb++){
            rhor[i] += (wfn[i + orb*bcSize]*wfn[i + orb*bcSize]);  
            //vrhor[i] = rhor[i]/femgrid[i+1];
        }
        rhor[i] = 2.0*rhor[i];
        rhor[i] = rhor[i]/femgrid[i+1];
        //printf("rho_r[%d] = %lf\n",i,rhor[i]);
    }
}
void ASFPM::getExchangeDens(double *rhox,int orb_a,int orb_b){
    for(int i=0; i<bcSize; i++){
        rhox[i] = 0.0;
        rhox[i] = wfn[orb_a*bcSize + i]*wfn[orb_b*bcSize + i];
        rhox[i] = rhox[i]/femgrid[i+1];
    }
}
void ASFPM::getDensityMatrix(double *densMat){ //May be implemented in the Atomic Structure Static Class
    int k=0;
    for(int i=0; i<bcSize; i++)
    {
        for(int j=0; j<bcSize; j++)
        {   densMat[k] = 0.0;
            for(int orb=0; orb<occOrb; orb++)
            {
                densMat[k] += 2.0*(wfn[i + orb*bcSize]*wfn[j + orb*bcSize]);
            }
            //densMat[k] = 2.0*densMat[k];
            k++;
        }
    }

}
// COULOMB TERM************************************
double * ASFPM::hartreePotential(){
    double *locdens = new double[bcSize];
    double *hpot = new double[bcSize];
    getDensity(locdens);
    solvePoissonEquation(hpot,locdens,numElectrons-1);
    delete []locdens;

    return hpot;
}
double * ASFPM::computeHartreePotential(){
    double *hpot = hartreePotential();

    double *vh_r = new double[total_nodes];
    for(int i=0; i<bcSize; i++){
        vh_r[i+1]= hpot[i]/femgrid[i+1];
    }
    vh_r[0]=doInterpolation(&femgrid[0],vh_r,2);
    vh_r[total_nodes-1]=(totQ-1)/femgrid[total_nodes-1];

    delete [] hpot;

    return vh_r;

}
// END COULOMB TERM************************************
// EXCHANGE TERM *************************************
double * ASFPM::exchangeAuxPotential(){
    double *aux = new double[bcSize];
    double *xdens = new double[bcSize];
    double *uxpot = new double[bcSize];
    for(int a=0; a<occOrb-1; a++){
        for(int b=1; b<occOrb; b++){
            if(a==b){
                continue;
            }
            getExchangeDens(xdens,a,b);
            solvePoissonEquation(aux,xdens,0.0);
            for(int i=0; i<bcSize; i++){
                uxpot[i] = 0.0;
                uxpot[i] = uxpot[i] + aux[i];
            }
        }
         
    }
    delete [] aux;
    delete [] xdens;

    return uxpot;
}
double * ASFPM::computeExchangePotential(){

    double *uxpot = exchangeAuxPotential();
    double *vh_r = new double[total_nodes];

    for(int i=0; i<bcSize; i++){
        vh_r[i+1]= uxpot[i]/femgrid[i+1];
    }
    vh_r[0]=doInterpolation(&femgrid[0],vh_r,2);
    vh_r[total_nodes-1]=0.0;

    delete [] uxpot;
    return vh_r;
}
//// END EXCHANGE TERM *************************************
double * ASFPM::computeHartreePotential(double *hpot,int rcIndex){
    
    double *vh_r = new double[total_nodes];

    for(int i=0; i<bcSize; i++){
        if(i>=rcIndex){
            vh_r[i] = 0.0;
        }
        else{
            vh_r[i+1]= hpot[i]/femgrid[i+1];
        }
    }
    vh_r[0]=doInterpolation(&femgrid[0],vh_r,2);
    vh_r[total_nodes-1]=0.0/femgrid[total_nodes-1];
    return vh_r;
}
double ASFPM::energyHF(double *hij, double *fij,double *densMat){
    
    double energy = 0.0;
    for(int i=0; i<bcSize*bcSize; i++)
    {
        energy = energy + 0.5*(densMat[i]*hij[i]);
        energy = energy + 0.5*(fij[i]*densMat[i]);
    }
    return 0.5*energy;

}
void ASFPM::singleDiagonalization(double *hij){
    std::cout<<"Single Diagonalization"<<std::endl;
    SumMatrices(&vij[0],&kij[0],hij,bcSize*bcSize);
    diag(bcSize,hij,&sij[0],eigenVal,wfn);
    //printf("Lowest Occupied Orbital Energy: %.10lf\n",0.5*eigenVal[0]);
}
void ASFPM::load(){
    //Create the fem grid:
    buildFemGrid(atomicN,r0,rN);
    if(atomicModel=="Soft-Walls"){
        int poiss_tot_nodes,pnodes;
        int rcIndex = femgrid.forcedInsertion(cutRad);
        femgrid.refineNear(rcIndex,0.001,8);
        poiss_tot_nodes = rcIndex+1;
        pnodes = poiss_tot_nodes-2;
        assambleMatricesFixedPoints(atomicN,rcIndex);
    }
    else{
        assambleMatricesFixedPoints(atomicN);
        //Hard-Walls and Free atomis
    }
    
}
void ASFPM::sfc(double *fij, double *hij, double *vhpot, double *vxpot,double &energy0){
    double *vhij = new double[bcSize*bcSize];
    double *vxij = new double[bcSize*bcSize];
    double *vtot = new double[bcSize*bcSize];
    double *densMat = new double[bcSize*bcSize];
    double *vXmat = new double[bcSize*bcSize];

    int iter = 0;
    double deltaEnergy = 100000000.0;
    double new_HF = energy0;
    double old_HF,orb_energy;
    printf("-----------------------------------------------------------\n");
    printf("-----------------------------------------------------------\n");
    printf("-------------- Self-Consistent Field ----------------------\n");
    printf("\n");
   // while(deltaEnergy>0.0000001){
    do{
        old_HF = new_HF;
        fixedPointsNumIntegration(vhij,vhpot);
        //fixedPointsNumIntegration(vxij,vxpot);
        //ScalarXMatrix(0.5,vxij,vXmat,bcSize,bcSize);
        //SumMatrices(vhij,vXmat,vtot,bcSize*bcSize);
        SumMatrices(hij,vhij,fij,bcSize*bcSize);
        diag(bcSize,fij,&sij[0],eigenVal,wfn);
        orb_energy = 0.5*eigenVal[0];
        wfnNormalization();
        getDensityMatrix(densMat);
        vhpot = computeHartreePotential();
        //vxpot = computeExchangePotential();
        new_HF = energyHF(hij,fij,densMat);
        printf("Hartree-Fock energy = %lf    orb = %lf   iteration: %d\n",new_HF,orb_energy,iter+1);
          
        deltaEnergy = fabs(new_HF-old_HF);
        iter++;
    }while(deltaEnergy>0.0000001);
    delete [] vhij;
    delete [] vxij;
    delete [] vtot;
    delete [] densMat;
    delete []vXmat;

}
//*** METHODS TO START THE PROGRAM *********

void ASFPM::startProgram(){
    //Checking initial data
    if(confType=="Soft-Walls"){
        if(atomicN!=1){
            double *uij{nullptr},*lij{nullptr};
            double *hij = new double[bcSize*bcSize];
            double *vhij = new double[bcSize*bcSize];
            double *densMat = new double[bcSize*bcSize];
            double *fij = new double[bcSize*bcSize];
            double *R_rho = new double[bcSize];
            double *R_hpot{nullptr};
            FillZeroMat(fij,bcSize,bcSize);
            FillZeroMat(vhij,bcSize,bcSize);
            FillZeroMat(hij,bcSize,bcSize);
            buildFemGrid(atomicN,r0,rN);
            int rcIndex = femgrid.forcedInsertion(cutRad);// Forces the Rc to appear in the grid
            femgrid.refineNear(rcIndex,0.001,8);
            //femgrid.refineNear(rcIndex,0.002,28);
             /* for(int i=0; i<total_nodes; i++){
                 printf("r[%d] = %lf\n",i,femgrid[i]);
            } */
            int k=0;
            for(int i=0; i<Ne; i++){
                for(int j=0; j<order+1; j++){
                    printf("E[%d].n[%d].x = %d = %lf delta\n",i,j,k, femgrid[k]);
                    k++; 
                }
                k--;
            }
            int poiss_tot_nodes = rcIndex+1;
            int pnodes = poiss_tot_nodes-2;
            //printf("Boundary Condition points for poisson: %d\n",pnodes);
            assambleMatricesFixedPoints(atomicN,rcIndex);// Must be modified 
            //assambleMatricesFixedPoints(atomicN);
            SumMatrices(&vij[0],&kij[0],hij,bcSize*bcSize);
            diag(bcSize,hij,&sij[0],eigenVal,wfn);
            printf("First orbital value: %.10lf\n",0.5*eigenVal[0]);
            wfnNormalization();
            getDensityMatrix(densMat);
            double energy0 = energyHF(hij,fij,densMat);
            printf("First Hartree-Fock Energy: %.10lf\n",energy0);
            getDensity();
            double qtot = computeTotalCharge(rcIndex);
            divideBy(R_rho,rho,&femgrid[0],bcSize);
            //printf("PINCHE PENDEJO REVISA ESTE VALOR DE LA CARGA TOTAL ----> %lf\n",totQ);
            printf("\n");
            printf("Solving the Poisson's equation to compute the Hartree Potential\n");
            double *hpot = new double[pnodes];
            double *bcVec  = new double[pnodes];
            //solvePoissonEquation(hpot,R_rho,totQ);
            lij = new double[pnodes*pnodes];
            uij = new double[pnodes*pnodes];
            
            assamblePoissonMatrices(lij,uij,bcVec,rcIndex);
            poissonSolver(hpot,uij,lij,R_rho,pnodes,bcVec,qtot);
            for(int i = 0; i<pnodes; i++){
                //printf("r[%d] = %lf   Rho(r) = %lf    UH(r) = %lf\n",i,femgrid[i+1],rho[i],hpot[i]);
                //printf("bvec = %lf\n",bcVec[i]);
            }
            R_hpot = computeHartreePotential(hpot,rcIndex);
              
            int iter = 0;
            double deltaEnergy = 100000000.0;
            double new_HF = energy0;
            double old_HF,orb_energy;
            while(deltaEnergy>0.0000001){
                old_HF = new_HF;
                fixedPointsNumIntegration(vhij,R_hpot);
                SumMatrices(hij,vhij,fij,bcSize*bcSize);
                diag(bcSize,fij,&sij[0],eigenVal,wfn);
                orb_energy = 0.5*eigenVal[0];
                wfnNormalization();
                getDensityMatrix(densMat);
                new_HF = energyHF(hij,fij,densMat);
                printf("Hartree-Fock energy = %lf    orb = %lf   iteration: %d\n",new_HF,orb_energy,iter+1);
                getDensity();
                qtot = computeTotalCharge(rcIndex);
                divideBy(R_rho,rho,&femgrid[0],bcSize);
                poissonSolver(hpot,uij,lij,R_rho,pnodes,bcVec,qtot);
                R_hpot = computeHartreePotential(hpot,rcIndex);
                deltaEnergy = fabs(new_HF-old_HF);
                iter++;
            }
            printf("-----------------------------------------------------------\n");
            printf("-----------------------------------------------------------\n");
            printf("--------------Hartree-Fock Final Results-------------------\n");
            printf("-----------------------------------------------------------\n");
            printf("Hatree-Fock Final Energy=%lf\n",new_HF);
            printf("-----------------------------------------------------------\n");
            printf("Occupied Orbital Energy = %lf\n",orb_energy);
            printf("-----------------------------------------------------------\n");
            printf("Final Charge = %lf\n",qtot);
            printf("-----------------------------------------------------------\n");
            for(int i=0; i<total_nodes; i++){
    
                printf("r[%d] = %lf   VH(r) = %lf\n",i,femgrid[i],R_hpot[i]);
            }     
            
            delete [] uij; 
            delete [] lij;
            delete [] hij;
            delete [] vhij;
            delete [] densMat;
            delete [] fij;
            delete [] R_rho;
            delete [] hpot;
            delete [] R_hpot;
            delete [] bcVec;

        }
        else{
            //We need to build the grid
            buildFemGrid(atomicN,r0,rN);
            //int pNodes = femgrid.forcedInsertion(cutRad);// Forces the Rc to appear in the grid
            assambleMatricesFixedPoints(atomicN);// Must be modified 
            double *hij = new double[bcSize*bcSize];
            int phase;
            int nodes = globalSize;
            SumMatrices(&vij[0],&kij[0],hij,bcSize*bcSize);
            diag(bcSize,hij,&sij[0],eigenVal,wfn);
            for(int i=0; i<total_nodes; i++){
                printf("r[%d] = %lf\n",i,femgrid[i]);
            } 
            printf("Highest Occupied Orbital Energy: %.10lf\n",0.5*eigenVal[0]);
            printf("Lowest Unoccupied Orbital Energy: %.10lf\n",0.5*eigenVal[1]);
            wfnNormalization();
            getDensity();
            //double qtot = computeTotalCharge();
            delete [] hij;
            }
    }
    else if(confType=="Hard-Walls"){
        printf("PUTO\n");
    }
    else{//FREE-ATOM
    }
}
void ASFPM::printWfn(int orbital){
    int realOrb = orbital-1;
    double wfnOrb;
    int phase;
    int phase2;
    std::ofstream wfnData;
    wfnData.open("wfn.dat",std::ios::out);
    if(!wfnData){
        std::cout<<"FILE NOT CREATED\n";
    } 
    getWfnPhase(fem_nodes,0,&phase,wfn);
    getWfnPhase(fem_nodes,1,&phase2,wfn);
    std::cout<<"The Wave Function Phase is: "<<phase<<std::endl;
    wfnData<<"r \t orb[0] \t orb[1] "<<std::endl;
    wfnData<<std::fixed<<std::setprecision(6)<<femgrid[0]<<"\t"<<0.0<<"\t"<<0.0<<std::endl;
    for(int i=0; i<fem_nodes; i++)
    {
        wfnData<<std::fixed<<std::setprecision(6)<<femgrid[i+1]<<"\t"<<wfn[i + fem_nodes*0]<<"\t"<<wfn[i + fem_nodes*1]<<std::endl;

    }
    //wfnData<<std::fixed<<std::setprecision(6)<<femgrid[0]<<"\t"<<0.0<<"\t"<<0.0<<std::endl;
    wfnData.close();
}














 /*    for(int i=0; i<total_nodes; i++){
    printf("Vr[%d] = %lf\n",i,vr[i]);
} */
/* for(int i=0; i<nodes-2; i++){
    for(int j=0; j<nodes-2; j++)
    {
        if(i==j)
        {
            printf("matV[%d] = %lf\n",i*(nodes-2) + j,vij[i*(nodes-2) + j]);
        }
    }
} */
//vij.printMatrix();
//femgrid.refineNear(pNodes,0.00001,2);
//femgrid.refineNear(pNodes,0.0002,3);
//femgrid.refineNear(pNodes,0.0002,3); //Works for order 2
//femgrid.printGrid(); //0.009 & 8 // 0.04 & 10) (0.009 & 3)
//sfemgrid.refineNear(pNodes,0.00,0);
//femgrid.printGrid();
//evaluateExternalPotential(vr,&femgrid[0],angular,atomicN,total_nodes,cutRad,wallValue);
//vr[0] = doInterpolation(&femgrid[0],vr,4);
//printf("Extrapolation Value v(r=0) = %lf\n",vr[0]);
/* int k=0;
for(int i=0; i<Ne; i++){
    for(int j=0; j<order+1; j++){
        printf("E[%d].n[%d].x = %d = %lf delta\n",i,j,k, femgrid[k]);
        k++; 
    }
    k--;
} */