#include "FemFixedPoints.hpp"

FixedPoints::FixedPoints(){

}
FixedPoints::FixedPoints(int _Ne, int _order, Grid<double> &inputGrid) 
    :Ne{_Ne},order{_order}{
    std::cout<<"Fixed Points FEM model Constructor\n";
    globalSize = Ne*order+1;
    poly=order + 1;
    bcSize = globalSize-2;
    linkmatSize = Ne*(order+1);
    linkMat.setMatrix(linkmatSize);
    buildLinkMAtrix();
    femgrid.setSize(globalSize);
    femgrid.copy(inputGrid);
    //bvec = new double[bcSize];
    //femgrid.printGrid();
}
FixedPoints::~FixedPoints(){
    delete [] lij;
    delete [] bvec;
};
void FixedPoints::buildLinkMAtrix(){
    for(int i=0; i<Ne; i++)
        {
            for(int j=0; j<poly; j++)
            {
                linkMat[poly*i+j] = poly*i+j-i;
            }
        }
}
void FixedPoints::extractBCvector(double *l_mat, double *boundaryVec,int nodes){
     for(int j=1; j<(nodes-1); j++)
        {
                boundaryVec[j-1] = l_mat[j + (j+1)*(nodes-1)];

        }

}
double *FixedPoints::getBoundaryConditionVector(){
    
    return bvec;   
}
void FixedPoints::reduceMatrix(double *matG, double *mat,int pts){
    int l=0;
    int n = globalSize; //Global Size for poisson problem
    for(int i=1;i<n-pts;i++)
    {
        for(int j=1;j<n-pts;j++)
        {
            mat[l] = matG[i*n+j]; //Overlap matrix
            l++;
        }
    }

}
void FixedPoints::setPoissonEquationData(double *matA,double *boundaryVec,int domain_size){
    lij = new double [domain_size*domain_size];
    bvec = new double[domain_size];
    poisson_domain_size = domain_size;
    for(int i=0; i<domain_size*domain_size; i++){
        lij[i] = matA[i];
    }
    for(int i=0; i<domain_size; i++){
        bvec[i] = boundaryVec[i];
    }
}
int FixedPoints::getBCSize(){
    return bcSize;
}
//Assamble the FEM Matrices
void FixedPoints::assambleMatrices(double *kij, double *sij, double *vij,double *boundaryVec,int atomicN){
    //FOR THE FIXED POINTS MODEL AND EXACT EXTERNAl POTENTIAL INTEGRATION
    double *s_matG = new double[globalSize*globalSize];
    double *k_matG = new double[globalSize*globalSize];
    double *v_matG = new double[globalSize*globalSize];
    FillZeroMat(v_matG,globalSize,globalSize);
    FillZeroMat(s_matG,globalSize,globalSize);
    FillZeroMat(k_matG,globalSize,globalSize);
    //**** ELEMENTAL MATRICES ***************
    double *feMatS{nullptr};
    double *feMatK{nullptr}; //This memory will be allocated by the function ;
    double *feMatV{nullptr};
    double x[poly]; //This method needs the points in the grid
    
    for(int ei=0; ei<Ne; ei++){
        //printf("e = %d\n",ei);
        for(int j=0; j<poly; j++){
            int indx = poly*ei + j; 
            int i = linkMat[indx];
            x[j] = femgrid[i];
            //printf("x[%d] = %lf\n",j,x[j]);
        }
        feMatS = getFixedPointsOverlapMatrices(x,order);
        feMatK = getFixedPointsKinectMatrices(x,order);
        feMatV = getFixedPointsAnaliticVr(x,order,atomicN);
        //if(x[])
        //printf("[0] = %lf\n",feMatV[0]);
        for(int nu=0; nu<poly; nu++){
            int index_nu = poly*ei + nu; 
            int l = linkMat[index_nu];
            for(int mu=0; mu<poly; mu++){

                int index_mu = poly*ei+mu;
                int m = linkMat[index_mu];
                s_matG[l*globalSize+m] += feMatS[poly*nu+mu];
                k_matG[l*globalSize+m] += feMatK[poly*nu+mu];
                //v_matG[l*globalSize+m] += feMatV[poly*nu+mu];
                //printf("v_matG[%d] = %lf + %lf = %lf\n",l*globalSize+m, v_matG[l*globalSize+m],feMatV[poly*nu+mu],v_matG[l*globalSize+m]);
                 
                    v_matG[l*globalSize+m] += feMatV[poly*nu+mu];
                    //printf("Vij[%d] = %lf\n",poly*nu+mu, feMatV[poly*nu+mu]);

                

            }
        }
    }
    extractBCvector(k_matG,boundaryVec,globalSize);
    
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
void FixedPoints::getSourceVector(double *sourceVec,double *matS, double *rightVec, double bcValue){
    double *aux = new double[poisson_domain_size];
    MatrixProduct(matS,rightVec,aux,poisson_domain_size,poisson_domain_size,1);
    for(int i = 0; i<poisson_domain_size; i++){
            sourceVec[i] = aux[i] - bvec[i]*bcValue;
            //printf("%lf\n",aux[i]);
    }
    delete [] aux;
}
//----- Poisson
void FixedPoints::solvePoissonEquation(double *xvector, double *sourceVec){
    int N = poisson_domain_size;
    int NRHS = 1;
    int LDA = N; 
    int LDB = N; 
    int ipiv[N];
    int info;
    

    double *right_vec = new double[poisson_domain_size];
    double *aux_vec = new double[poisson_domain_size];
    double *aux_mat = new double[poisson_domain_size*poisson_domain_size];  

    ScalarXMatrix(1.0,lij,aux_mat,poisson_domain_size,poisson_domain_size);

    //MatrixProduct(&sij[0],rho_r,right_vec,bcSize,bcSize,NRHS);
    for(int i=0; i<bcSize; i++){
        aux_vec[i] = sourceVec[i];
        //qtot  = qtot + right_vec[i];
        //printf("right_vec = %lf - %lf\n",right_vec[i],-bvec[i]);
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
        xvector[i] = aux_vec[i];
            //printf("HP[%d] = % lf\n",i,aux_vec[i]);
    }

    delete [] aux_mat;
    delete [] right_vec;
    delete [] aux_vec;
}
void FixedPoints::femNumericalIntegration(double *inMat,double *pot){

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
    // it eliminates 1 point of the domain

    //**********************************************************

    // FREE ALL THE SHIT ALLOCATED ABOVE
    delete [] feMatP;
    delete [] p_matG;

}
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