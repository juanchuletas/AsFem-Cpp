#if !defined(_EXTERNAl_POTENTIAL_H_)
#define _EXTERNAl_POTENTIAL_H_

#include<iostream>
#include<cmath>

void evaluateExternalPotential(double *vr,double *x,int angular,int atomicN,int points);
void evaluateExternalPotential(double *vr,double *x,int angular,int atomicN,int Ne,double cutRad, double wallValue);
void evaluateExternalPotential(double *vr,double *x,int angular,int atomicN,int Ne,double screen);

#endif // _EXTERNAl_POTENTIAL_H_


