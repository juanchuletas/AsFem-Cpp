#include "FiniteElement.hpp"




FEM::FEM(){
    Ne = 0; 
    order = 0;
    globalSize = 0;
    poly = 0; 
    bcSize = 0; 
    linkmatSize = 0; 
}
FEM::FEM(int _Ne, int _order, std::string _femModel, std::string _gridName)
:Ne{_Ne},order{_order},femModel{_femModel},gridType{_gridName}{
    std::cout<<"Constructor 1. \n";
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
FEM::FEM(int _Ne, int _order, int _poissonNe, std::string _femModel,std::string gridName)
 :Ne{_Ne},order{_order},femModel{_femModel},gridType{gridName}{
    poly = order + 1;
    globalSize = Ne*order+1;
    bcSize = globalSize-2; //Dirichet Boundary Conditions
    poissNodes = (poissNodes*order+1);
    linkmatSize = Ne*(order+1);
    linkMat.setMatrix(linkmatSize);
    buildLinkMAtrix();
    sij.setMatrix(bcSize*bcSize);
    sij.fillWithZeros();
    vij.setMatrix(bcSize*bcSize);
    vij.fillWithZeros();
    kij.setMatrix(bcSize*bcSize);
    kij.fillWithZeros();
    //femgrid.setGridData(0.0, 20.0, 10, 2,"Froese-Fischer")
}
FEM::~FEM(){
    std::cout<<"FEM destructor called\n";
    delete [] bvec;
}
void FEM::buildLinkMAtrix(){
	for(int i=0; i<Ne; i++)
        {
            for(int j=0; j<poly; j++)
            {
                linkMat[poly*i+j] = poly*i+j-i;
            }
        }
}
void FEM::setFemData(int inNe,int inOrder){
    order = inOrder;
    Ne = inNe;
    poly = order+1;
    linkmatSize = Ne*(order+1);
    linkMat.setMatrix(linkmatSize);
    buildLinkMAtrix();
}
int & FEM::getLinkMatIndex(int i){
 return linkMat[i];
}
void FEM::buildFemGrid(int atomicN,double r0,double rN){
    if(femModel=="Fixed Elements"){
        femgrid.setGridData(r0, rN, Ne, order,gridType, atomicN);
        femgrid.createGrid();
    }
    else{//FIXED POINTS MODEL
        femgrid.setGridData(r0, rN, Ne, order,gridType, atomicN);
        femgrid.createGrid(femModel);
    }

}
//PRIVATE METHOD***
void FEM::extractBCvector(double *l_mat,int nodes){
     for(int j=1; j<(nodes-1); j++)
        {
                bvec[j-1] = l_mat[j + (j+1)*(nodes-1)];

        }

}
void FEM::reduceMatrix(double *matG, double *mat,int size){
    int l=0;
    int n = size; //Global Size for poisson problem
    for(int i=1;i<n-1;i++)
    {
        for(int j=1;j<n-1;j++)
        {
            mat[l] = matG[i*n+j]; //Overlap matrix
            l++;
        }
    }

}
void FEM::reduceMatrix(double *matG, Matrix<double> &mat,int size){
    int l=0;
    int n = size; //Global Size for poisson problem
    for(int i=1;i<n-1;i++)
    {
        for(int j=1;j<n-1;j++)
        {
            mat[l] = matG[i*n+j]; //Overlap matrix
            l++;
        }
    }

}
void FEM::assambleMatricesFixedPoints(int atomicN){
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
            }
        }
    }
    reduceMatrix(s_matG,sij,globalSize);
    reduceMatrix(k_matG,kij,globalSize);
    reduceMatrix(v_matG,vij,globalSize);
    
    //**********************************************************

    // FREE ALL THE SHIT ALLOCATED ABOVE
    delete [] feMatS;
    delete [] feMatK;
    delete [] feMatV;
    delete [] s_matG;
    delete [] k_matG;
    delete [] v_matG;
    
}

//**** POISON SOLVER *****************************************
void FEM::solvePoissonEquation(double *hpot, double *rho_r,double hp){
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
        //printf("f = %lf\n",rho[i]);
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

void FEM::assambleMatricesFixedElements(double *pot){
    //FOR THE FIXED ELEMENTS MODEL AND NUMERIC INTEGRATION OF THE EXTERNAL POTENTIAL
    //****** GLOBAL MATRICES *************
    double *s_matG = new double[globalSize*globalSize];
    double *k_matG = new double[globalSize*globalSize];
    double *v_matG = new double[globalSize*globalSize];
    //**** ELEMENTAL MATRICES ***************
    double *feMatS;
    double *feMatK; //This memory will be allocated by the function ;
    double *feMatV;
    double v[poly];
    //**************************************

    std::cout<<femModel<<std::endl;
    feMatS = ElementalOverlapMatrix(order);
    feMatK = ElementalKinectMatrix(order);
    double coeff;
    for(int ei=0; ei<Ne; ei++)
    {
        coeff = 0.5*femgrid.getElementSize(ei);
        for(int j=0; j<poly; j++){
            int indx = poly*ei+j;
            int i  = getLinkMatIndex(indx);
            v[j] = pot[i];
        }
        feMatV = ElementalPotentialMatrix(order,v);
        for(int nu=0; nu<poly; nu++)
        {
            int index_nu = poly*ei + nu;
            int l = getLinkMatIndex(index_nu);
            for(int mu = 0; mu<poly; mu++)
            {
                int index_mu = poly*ei+mu;
                int m = getLinkMatIndex(index_mu);
                s_matG[l*globalSize+m] += coeff*feMatS[poly*nu+mu];
                v_matG[l*globalSize+m] += coeff*feMatV[poly*nu+mu];
                k_matG[l*globalSize+m] += feMatK[poly*nu+mu]/coeff;
            }
        }
    }
    extractBCvector(k_matG,globalSize);
    //**** This part reduce the global matrices according to the 
    // proper boundary conditions***********************
    reduceMatrix(s_matG,sij,globalSize);
    reduceMatrix(k_matG,kij,globalSize);
    reduceMatrix(v_matG,vij,globalSize);

    //**********************************************************

    // FREE ALL THE SHIT ALLOCATED ABOVE
    delete [] feMatS;
    delete [] feMatK;
    delete [] feMatV;
    delete [] s_matG;
    delete [] k_matG;
    delete [] v_matG;
}
//************* METHODS FOR NUMERICAL INTEGRATION
void FEM::fixedPointsNumIntegration(double *inMat,double *pot){

    double *p_matG = new double[globalSize*globalSize];
    FillZeroMat(p_matG,globalSize,globalSize);
    double *feMatP;
    double v[poly];
    double x[poly];
    
    for(int ei=0; ei<Ne; ei++){
        for(int j=0; j<poly; j++){
            int indx = poly*ei+j;
            int i = getLinkMatIndex(indx);
            v[j] = pot[i];
            x[j] = femgrid[i];
        }
        feMatP = getFixedPointsNumPotMatrices(x,v,order);
        for(int nu=0; nu<poly; nu++) {
            int index_nu = poly*ei + nu;
            int l = getLinkMatIndex(index_nu);
            for(int mu = 0; mu<poly; mu++){
                int index_mu = poly*ei+mu;
                int m = getLinkMatIndex(index_mu);
                p_matG[l*globalSize+m] += feMatP[poly*nu+mu];
                //printf("fem[%d] = %lf \n",l*globalSize+m,p_matG[l*globalSize+m]);
            }
        }
    }
    //**** This part reduce the global matrices according to the 
    // proper boundary conditions***********************
    reduceMatrix(p_matG,inMat,globalSize);

    //**********************************************************

    // FREE ALL THE SHIT ALLOCATED ABOVE
    delete [] feMatP;
    delete [] p_matG;

}
void FEM::fixedElementsNumIntegration(double *inMat, double *pot)
{
    double *p_matG = new double[globalSize*globalSize];
    FillZeroMat(p_matG,globalSize,globalSize);
    double *feMatP;
    double v[poly];
    double coeff;
    for(int ei=0; ei<Ne; ei++)
    {
        coeff = 0.5*femgrid.getElementSize(ei);
        for(int j=0; j<poly; j++){
            int indx = poly*ei+j;
            int i  = getLinkMatIndex(indx);
            v[j] = pot[i];
            //printf("v[%d] = %lf\n",j,v[j]);
        }
        feMatP = ElementalPotentialMatrix(order,v);
        for(int nu=0; nu<poly; nu++)
        {
            int index_nu = poly*ei + nu;
            int l = getLinkMatIndex(index_nu);
            for(int mu = 0; mu<poly; mu++)
            {
                int index_mu = poly*ei+mu;
                int m = getLinkMatIndex(index_mu);
                p_matG[l*globalSize+m] += coeff*feMatP[poly*nu+mu];
                //printf("fem[%d] = %lf \n",l*globalSize+m,p_matG[l*globalSize+m]);
            }
        }
    }
    //**** This part reduce the global matrices according to the 
    // proper boundary conditions***********************
    reduceMatrix(p_matG,inMat,globalSize);

    //**********************************************************

    // FREE ALL THE SHIT ALLOCATED ABOVE
    delete [] feMatP;
    delete [] p_matG;
 
}
void FEM::assamblePoissonMatrices(double *lij, double *uij,int pNe){
    int pNodes = pNe*order+1; //pNodes->poisson Nodes
    int pbcNodes = pNodes-2;
    double *poiss_matL = new double[pNodes*pNodes];
    double *poiss_matU= new double[pNodes*pNodes];
    FillZeroMat(poiss_matL,pNodes,pNodes);
    FillZeroMat(poiss_matU,pNodes,pNodes);
    double *feMatL, *feMatU;

    double x[poly]; //This method needs the points in the grid
    for(int ei=0; ei<pNe; ei++){
        for(int j=0; j<poly; j++){
            int indx = poly*ei + j; 
            int i = linkMat[indx];
            x[j] = femgrid[i];
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
                poiss_matU[l*pNodes+m] += feMatU[poly*nu+mu];
                poiss_matL[l*pNodes+m] += feMatL[poly*nu+mu];
            }
        }
    }
    reduceMatrix(poiss_matU,uij,pNodes);
    reduceMatrix(poiss_matL,lij,pNodes);
    delete [] poiss_matL;
    delete [] poiss_matU;

}

// *** END METHODS
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