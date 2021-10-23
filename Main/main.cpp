#include<iostream>
//#include "../ElectronicStructure/closedShell.hpp"
//#include "../Atomic_Computations/atomicStructure.hpp"
#include "../ElectronicStructure/closedShell.hpp"

int main(){
	grid_tools::froese_fischer::hmf = 1.0/32.0;
	grid_tools::froese_fischer::rmf = 5.0;
	int Ne;
	int order = 3;
	double r0 = 0;
	double rc = 2.0;
	double rinfty = 250.0;
	int atomicN = 4; 
	int charge = 0;
	int virtualOrbs = 1;
	int numElec = atomicN-charge; 
	double totpoints = grid_tools::froese_fischer::inverseKernel(rinfty,atomicN);
	int points = floor(totpoints);
	std::cout<<"Points gained: "<<points<<std::endl;
	if(points%2==0){
		printf("True\n");
		points = points - 1;
	}
	//points = 100;
	Ne = (points-1)/order;
	int pts = Ne*order+1;
	std::cout<<"Elements for this run: "<<Ne<<std::endl;
	std::cout<<"Points for this run: "<<Ne*order+1<<std::endl;
	std::cout<<"Order: "<<order<<std::endl;
	std::string mesh = "Froese-Fischer";
	std::string mesh2 = "Chebyshev";
	std::string kind = "Fixed-Points";
	//Grid<double> myGrid;
	Grid<double> myGrid{r0,rinfty,Ne,order,"Exponential"};
	myGrid.froeseFischer(atomicN);
	//myGrid.chebyshev();
	//myGrid.printGrid();
	//myGrid.setGridData(r0,rinfty,Ne,order,mesh,atomicN);
	//myGrid.createGrid(kind);
	//printf("SIZE 1: %d\n",myGrid.size());
	ClosedShell::ElectronicStructure asfem{Ne,order,myGrid,atomicN,numElec,virtualOrbs,rinfty};
/* 	double *matS = new double[(order+1)*(order+1)];
	double x[] ={myGrid[0],myGrid[1],myGrid[2],myGrid[3],myGrid[4],myGrid[5],myGrid[6]}; 
	matS = getFixedPointsOverlapMatrices(x,order);
	double x0,x1,x2,x3;
	x0 = x[0]; x1=x[1]; x2=x[2]; x3=x[3];
	for(int i = 0; i<(6+1); i++){
		printf("x[%d] = %.10lf\n",i,x[i]);
	}
	double  res = (-1.0/210.0)*(pow(x0 - x3,5)*(2*pow(x0,2) - 7*x0*x2 + 7*pow(x2,2) + 3*x0*x3 - 7*x2*x3 + 2*pow(x3,2)))/(pow(x0 - x1,2)*pow(x1 - x2,2)*pow(x1 - x3,2));
	printf("res = %lf\n",res);
	for(int i=0; i<(order+1)*(order+1); i++){
		
	} */
	int M,N;
	M = asfem.getBCSize();
	N = asfem.getBCSize();
	int globalSize = Ne*order + 1; 

	double *k_mat = new double[M*N];
	double *s_mat = new double[M*N];
	double *hcore_mat = new double[M*N];
	double *vne_mat = new double[M*N];
	double *matCoeffs = new double[M*N];
	double *eigenVec = new double[M];

	//---- Boundary condition vector
	//---  this vector is necessary for the poisson equation 
	double *bcVector = new double[M];
	//------------------------------------------------------
	// Free-Model FEM matrices *****
	asfem.assambleMatrices(k_mat,s_mat,vne_mat,bcVector,atomicN);
	asfem.setPoissonEquationData(k_mat,bcVector,asfem.getBCSize());
	bool fullSolver = false;
	int tol = 3;
	SumMatrices(k_mat,vne_mat,hcore_mat,M*N);	
	//asfem.eigenSystemSCF(hcore_mat,s_mat,matCoeffs,eigenVec,fullSolver, tol); 
	asfem.iterativeSCF(hcore_mat,s_mat,matCoeffs);	
	delete [] k_mat;
	delete [] s_mat;
	delete [] vne_mat;
	delete [] matCoeffs;
	delete [] eigenVec;  
	return 0;
}
