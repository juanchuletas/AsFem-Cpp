#include "FiniteElementFixedPoints.hpp"


FEMFP::FEMFP(){
    Ne = 0; 
    order = 0;
    globalSize = 0;
    poly = 0; 
    bcSize = 0; 
    linkmatSize = 0; 
}
FEMFP::FEMFP(int _Ne, int _order, std::string _femModel, std::string _gridName)
:Ne{_Ne},order{_order},femModel{_femModel},gridType{_gridName}{
    //This constructor works when you have to solve the Poisson Equation with the same
    // Size of the eigenvalue size
    poly=order + 1;
    globalSize = Ne*order + 1;
    bcSize = globalSize-2;
    linkmatSize = Ne*(order+1);
    linkMat.setMatrix(linkmatSize);
    buildLinkMAtrix();
    sij.setMatrix(bcSize*bcSize);
    sij.fillWithZeros();
    vij.setMatrix(bcSize*bcSize);
    vij.fillWithZeros();
    kij.setMatrix(bcSize*bcSize);
    kij.fillWithZeros();
    bvec = new double[bcSize];
    
}
FEMFP::~FEMFP(){
    std::cout<<"FEM destructor called\n";
    delete [] bvec;
}
//***** Methods ***********

void FEMFP::buildLinkMAtrix(){
    for(int i=0; i<Ne; i++)
        {
            for(int j=0; j<poly; j++)
            {
                linkMat[poly*i+j] = poly*i+j-i;
            }
        }
}
void FEMFP::buildFemGrid(int atomicN,double r0,double rN){

    femgrid.setGridData(r0, rN, Ne, order,gridType, atomicN);
    femgrid.createGrid(femModel);
}
void FEMFP::assambleMatricesFixedPoints(int atomicN){
    //FOR THE FIXED POINTS MODEL AND EXACT EXTERNAl POTENTIAL INTEGRATION
    double *s_matG = new double[globalSize*globalSize];
    double *k_matG = new double[globalSize*globalSize];
    double *v_matG = new double[globalSize*globalSize];
    //**** ELEMENTAL MATRICES ***************
    double *feMatS;
    double *feMatK; //This memory will be allocated by the function ;
    double *feMatV;
    double x[poly]; //This method needs the points in the grid
    for(int ei=0; ei<Ne; ei++){
        for(int j=0; j<poly; j++){
            int indx = poly*ei + j; 
            int i = linkMat[indx];
            x[j] = femgrid[i];
        }
        feMatS = getFixedPointsOverlapMatrices(x,order);
        feMatK = getFixedPointsKinectMatrices(x,order);
        feMatV = getFixedPointsAnaliticVr(x,order,atomicN);
        //printf("[0] = %lf\n",feMatV[0]);
        for(int nu=0; nu<poly; nu++){
            int index_nu = poly*ei + nu; 
            int l = linkMat[index_nu];
            for(int mu=0; mu<poly; mu++){

                int index_mu = poly*ei+mu;
                int m = linkMat[index_mu];
                s_matG[l*globalSize+m] += feMatS[poly*nu+mu];
                k_matG[l*globalSize+m] += feMatK[poly*nu+mu];
                v_matG[l*globalSize+m] += feMatV[poly*nu+mu];
                //printf("s_matG = %lf\n",s_matG[l*globalSize+m]);
            }
        }
    }
    extractBCvector(k_matG,globalSize);
    
    reduceMatrix(s_matG,sij,1);
    reduceMatrix(k_matG,kij,1);
    reduceMatrix(v_matG,vij,1);
    
    //**********************************************************

    // FREE ALL THE SHIT ALLOCATED ABOVE
    delete [] feMatS;
    delete [] feMatK;
    delete [] feMatV;
    delete [] s_matG;
    delete [] k_matG;
    delete [] v_matG;
    
}
// For numerical integrations***********
void FEMFP::fixedPointsNumIntegration(double *inMat,double *pot){

    double *p_matG = new double[globalSize*globalSize];
    FillZeroMat(p_matG,globalSize,globalSize);
    double *feMatP;
    double v[poly];
    double x[poly];
    for(int ei=0; ei<Ne; ei++){
        for(int j=0; j<poly; j++){
            int indx = poly*ei+j;
            int i = linkMat[indx];
            v[j] = pot[i];
            x[j] = femgrid[i];
        }
        feMatP = getFixedPointsNumPotMatrices(x,v,order);
        for(int nu=0; nu<poly; nu++) {
            int index_nu = poly*ei + nu;
            int l = linkMat[index_nu];
            for(int mu = 0; mu<poly; mu++){
                int index_mu = poly*ei+mu;
                int m = linkMat[index_mu];
                p_matG[l*globalSize+m] += feMatP[poly*nu+mu];
                //printf("fem[%d] = %lf \n",l*globalSize+m,p_matG[l*globalSize+m]);
            }
        }
    }
    //**** This part reduce the global matrices according to the 
    // proper boundary conditions***********************
    reduceMatrix(p_matG,inMat,1);

    //**********************************************************

    // FREE ALL THE SHIT ALLOCATED ABOVE
    delete [] feMatP;
    delete [] p_matG;

}
// POISSON EQUATION SOLVER***************
void FEMFP::solvePoissonEquation(double *hpot, double *rho_r,double hp){
    int N = bcSize;
    int NRHS = 1;
    int LDA = N; 
    int LDB = N; 
    int ipiv[N];
    int info;
    

    double *right_vec = new double[bcSize];
    double *aux_vec = new double[bcSize];
    double *aux_mat = new double[bcSize*bcSize];  

    ScalarXMatrix(1.0,&kij[0],aux_mat,bcSize,bcSize);

    MatrixProduct(&sij[0],rho_r,right_vec,bcSize,bcSize,NRHS);
    for(int i=0; i<bcSize; i++){
        aux_vec[i] = right_vec[i]  - bvec[i]*hp;
        //qtot  = qtot + right_vec[i];
        printf("right_vec = %lf - %lf\n",right_vec[i],-bvec[i]);
    }
    dgesv_(&N,&NRHS,aux_mat,&LDA,ipiv,aux_vec,&LDB,&info);
    if( info > 0 ) 
    {
            printf( "The diagonal element of the triangular factor of A,\n" );
            printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
            printf( "the solution could not be computed.\n" );
            exit(1);
    }
    for(int i=0; i<bcSize; i++)
    {
        hpot[i] = aux_vec[i];
            //printf("HP[%d] = % lf\n",i,aux_vec[i]);
    }

    delete [] aux_mat;
    delete [] right_vec;
    delete [] aux_vec;
}
void FEMFP::assamblePoissonMatrices(double *lij, double *uij,double *bcVec,int rcIndex){
    int nele = (rcIndex-1)/order;
    nele = nele + 1;
    int poiss_points = rcIndex+1;
    int poiss_bc = poiss_points-2;
    double *poiss_matL = new double[poiss_points*poiss_points];
    double *poiss_matU= new double[poiss_points*poiss_points];
    FillZeroMat(poiss_matL,poiss_points,poiss_points);
    FillZeroMat(poiss_matU,poiss_points,poiss_points);
    double *feMatL, *feMatU;
    double x[poly]; //This method needs the points in the grid
    for(int ei=0; ei<nele; ei++){
        //printf("e = %d\n",ei);
        for(int j=0; j<poly; j++){
            int indx = poly*ei + j; 
            int i = linkMat[indx];
            x[j] = femgrid[i];
            //printf("x[%d] = %lf\n",j,x[j]);
        }
        feMatU = getFixedPointsOverlapMatrices(x,order);
        feMatL = getFixedPointsKinectMatrices(x,order);
        //printf("[0] = %lf\n",feMatV[0]);
        for(int nu=0; nu<poly; nu++){
            int index_nu = poly*ei + nu; 
            int l = linkMat[index_nu];
            for(int mu=0; mu<poly; mu++){

                int index_mu = poly*ei+mu;
                int m = linkMat[index_mu];
                poiss_matU[l*poiss_points+m] += feMatU[poly*nu+mu];
                poiss_matL[l*poiss_points+m] += feMatL[poly*nu+mu];
                //printf("Lij[%d] = %lf\n",poly*nu+mu, poiss_matL[l*globalSize+m]);
            }
        }
    }
    extractBCvector(poiss_matL,bcVec,poiss_points);
    reduceMatrix(poiss_matU,uij,poiss_points,1);
    reduceMatrix(poiss_matL,lij,poiss_points,1);
    //printf("matU[0] = %lf\n",uij[0]);
    delete [] poiss_matL;
    delete [] poiss_matU;

}
// END METHODS
double integrateElement(int ei,int order,double *feMatS, int *link_mat, double coeff, double *cf){
    int poly = order +1; 
    double phi[poly];
    for(int j=0; j<poly; j++){
        phi[j] = cf[link_mat[poly*ei+j]];
    }
    double value = 0.0;
    for(int mu=0; mu<poly; mu++){
        for(int nu=0; nu<poly; nu++){
            value = value + phi[mu]*phi[nu]*(coeff*feMatS[poly*mu + nu]);
        }
    }
    return value;

}