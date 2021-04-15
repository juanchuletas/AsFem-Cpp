#if !defined(_MATRICES_OPERATIONS)
#define _MATRICES_OPERATIONS
#include<iostream>
#include<cmath>
void SumMatrices(double *A, double *B, double *C, int size);
void MatrixProduct(double *mat_A,double *mat_B, double *mat_C,int N,int M, int P);
void MatXvec(double *matA,double *vec,double *res,int SIZE);
void ScalarXMatrix(double coeff,double *mat_A,double *mat_Res,int N,int M);
void ColumnMayor(double *mat_A, double *mat_C, int N, int M);
void RowMayor(double *mat_A, double *mat_C, int N, int M);
#endif // _MATRICES_OPERATIONS
