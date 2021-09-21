#if !defined(_MATH_METHODS_H_)
#define _MATH_METHODS_H_
#include<iostream>
#include<cmath>
#include "matrices_operations.hpp"
extern "C" void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
extern "C" void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
void poissonSolver(double *hartree_vec,double *sij,double *lij,double *rho,int pnodes,double *b,double hp);
void diag (int n, double *h, double *s, double *e, double *v);
double doInterpolation(double *xin,double *yin,int grade);
double nPolExtrapolation(double *x,double *y, int N,double target);
void getWfnPhase(int nodes, int orb,int *phase,double *wfn);
#endif // _MATH_METHODS_H_
