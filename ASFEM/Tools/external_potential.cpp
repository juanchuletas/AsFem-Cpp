#include "external_potential.hpp"
void evaluateExternalPotential(double *vr,double *x,int angular,int atomicN,int points){
    double z = static_cast<double> (atomicN);
    double l = static_cast<double>(angular);
    vr[0] = 0.f;
    for(int i=0; i<points; i++){
        double target  = ((l*(l+1))/x[i]*x[i] - 2.0*z/x[i]);
        vr[i] = target;
    }
}
void evaluateExternalPotential(double *vr,double *x,int angular,int atomicN,int Ne,double cutRad, double wallValue){
    double z = static_cast<double> (atomicN);
    double l = static_cast<double>(angular);
    for(int i=0; i<Ne; i++){
        if(x[i]<cutRad){
            double target  = ((l*(l+1))/x[i]*x[i] - 2.0*z/x[i]);
            vr[i] = target;
        }
        else{
            vr[i] = wallValue;
        }
    }
}
void evaluateExternalPotential(double *vr,double *x,int angular,int atomicN,int Ne,double screen){
    double z = static_cast<double> (atomicN);
    double l = static_cast<double>(angular);
    double b = screen;
    for(int i=0; i<Ne; i++){
        double target = ((l*(l+1))/x[i]*x[i] - ((2.0*z*exp(-x[i]*b))/x[i]));
        vr[i] = target;
    }

}
