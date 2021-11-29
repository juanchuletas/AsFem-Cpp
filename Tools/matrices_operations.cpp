#include "matrices_operations.hpp"


void SumMatrices(double *A, double *B, double *C, int size)
{
	for(int i=0; i<size; i++){
	
		C[i] = A[i]+B[i];
	}

}

void FillZeroMat(double *mat,int M,int N)
{
	for(int i=0; i<M*N; i++)
	{
		mat[i] = 0.0;
	}

}

void MatrixProduct(double *mat_A,double *mat_B, double *mat_C,int N,int M, int P)
{

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<P; j++)
		{
			mat_C[i*P + j] = 0.0;
			for(int k=0; k<M; k++)
			{
				mat_C[i*P + j] += mat_A[i*M + k]*mat_B[k + j*P];
				
			}
		}
	} 

}
void matMult(double *A, int rowsA,int colsA, double *B,int rowsB, int colsB,double *C){

  double sum;
  for(int i=0; i<rowsA; i++){
	for(int j=0; j<colsB; j++){
	  C[i*colsA+j] = 0.0;
	  //printf("C[%d]\n",i*colsA+j);
	  for(int k=0; k<rowsB; k++){
		//printf("sum[%d]= A[%d]*B[%d]\n",i*colsA+j,  i*colsA + k,   j*colsB+k);
	    C[i*colsA+j] += A[i*colsA + k]*B[j*colsB+k];
	  }
      
	}
  }
}
void MatXvec(double *matA,double *vec,double *res,int SIZE)
{
	for(int i=0; i<SIZE; i++)
	{
		for(int j=0; j<SIZE; j++)
		{
			res[i] += matA[j + i*SIZE]*vec[j];
		}
	}
}
void ScalarXMatrix(double coeff,double *mat_A,double *mat_Res,int N,int M)
{
	for(int i=0; i<N*M; i++)
	{
		mat_Res[i] = coeff*mat_A[i];
	}
}
void ColumnMayor(double *mat_A, double *mat_C, int N, int M)
{
	//printf("MatC:\tMatA\n");
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<M; j++)
		{
			mat_C[j*M + i] = mat_A[i*N + j];
			 
		}
	}
}
void RowMayor(double *mat_A, double *mat_C, int N, int M)
{
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<M; j++)
		{
			mat_C[i*N + j] = mat_A[j*M + i];
			 
		}
	}
}
void divideBy(double *newrho, double *rho,double *femnodes,int nfembasis){


	double factor;
	int l=0;
	for(int i=0; i<nfembasis; i++){
		
		newrho[i] = rho[i]/femnodes[i+1];
	
	}

}
