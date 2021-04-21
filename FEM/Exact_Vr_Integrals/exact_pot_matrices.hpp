#include<iostream>
#include<cmath>
#if !defined(_EXACT_Vr_INTEGRALS_H_)
#define _EXACT_Vr_INTEGRALS_H_

double *getExactPotMatrices(double *x,int order,int atomicN)
{
	int P = order +1;
	double *matVne = new double[P*P];
	double z = static_cast<double>(atomicN);
	switch(order)
	{
		case 2:
		{
			double x0,x1,x2;
			x0=x[0]; x1=x[1]; x2=x[2];
			if(fabs(x0)==0.f){
				matVne[0] = 0.0;
				matVne[1] = 0.0;
				matVne[2] = 0.0;
				matVne[3] = 0.0;
				matVne[4] = z*(-0.16666666666666666*pow(x2,4)/(pow(x1,2)*pow(x1 - x2,2)));
				matVne[5] = z*(((2*x1 - x2)*pow(x2,2))/(6.*x1*pow(x1 - x2,2)));
				matVne[6] = 0.0;
				matVne[7] = z*(((2*x1 - x2)*pow(x2,2))/(6.*x1*pow(x1 - x2,2)));
				matVne[8] = z*((-6*pow(x1,2) + 8*x1*x2 - 3*pow(x2,2))/(6.*pow(x1 - x2,2)));

			}
			else{
				matVne[0] = z*(((x0 - x2)*(3*pow(x0,3) - pow(x0,2)*(8*x1 + 5*x2) + x2*(-18*pow(x1,2) - 8*x1*x2 + pow(x2,2)) + x0*(6*pow(x1,2) + 16*x1*x2 + pow(x2,2))) + 12*pow(x1,2)*pow(x2,2)*log(x0) - 12*pow(x1,2)*pow(x2,2)*log(x2))/(6.*pow(x0 - x1,2)*pow(x0 - x2,2)));
				matVne[1] = z*(((x0 - x2)*(pow(x0,3) + (4*x1 - x2)*pow(x2,2) - pow(x0,2)*(2*x1 + 3*x2) + x0*x2*(10*x1 + 3*x2)) - 12*x0*x1*pow(x2,2)*log(x0) + 12*x0*x1*pow(x2,2)*log(x2))/(6.*pow(x0 - x1,2)*(x0 - x2)*(x1 - x2)));
				matVne[2] = z*((-1.0/6.0)*((x0 - x2)*(pow(x0,3) - pow(x0,2)*(4*x1 + x2) + x0*(6*pow(x1,2) + 8*x1*x2 - pow(x2,2)) + x2*(6*pow(x1,2) - 4*x1*x2 + pow(x2,2))) - 12*x0*pow(x1,2)*x2*log(x0) + 12*x0*pow(x1,2)*x2*log(x2))/((x0 - x1)*pow(x0 - x2,2)*(x1 - x2)));
				matVne[3] = matVne[1];
				matVne[4] = z*((pow(x0,4) - 8*pow(x0,3)*x2 + 8*x0*pow(x2,3) - pow(x2,4) + 12*pow(x0,2)*pow(x2,2)*log(x0) - 12*pow(x0,2)*pow(x2,2)*log(x2))/(6.*pow(x0 - x1,2)*pow(x1 - x2,2)));
				matVne[5] = z*((-1.0/6.0)*((x0 - x2)*(pow(x0,3) + (2*x1 - x2)*pow(x2,2) + x0*x2*(-10*x1 + 3*x2) - pow(x0,2)*(4*x1 + 3*x2)) + 12*pow(x0,2)*x1*x2*log(x0) - 12*pow(x0,2)*x1*x2*log(x2))/((x0 - x1)*(x0 - x2)*pow(x1 - x2,2)));
				matVne[6] = matVne[2];
				matVne[7] = matVne[5];
				matVne[8] = z*(((x0 - x2)*(pow(x0,3) + pow(x0,2)*(-8*x1 + x2) + x0*(-18*pow(x1,2) + 16*x1*x2 - 5*pow(x2,2)) + x2*(6*pow(x1,2) - 8*x1*x2 + 3*pow(x2,2))) + 12*pow(x0,2)*pow(x1,2)*log(x0) - 12*pow(x0,2)*pow(x1,2)*log(x2))/(6.*pow(x0 - x2,2)*pow(x1 - x2,2)));

			}

			
			break;
		}
		case 3:
		{
			double x0,x1,x2,x3;
			x0 = x[0]; x1=x[1]; x2=x[2]; x3=x[3];
			if(x0==0.f){
				matVne[0] = 0.0;
				matVne[1] = 0.0;
				matVne[2] = 0.0;
				matVne[3] = 0.0;
				matVne[4] = 0.0;
				matVne[5] = matVne[5] = -0.03333333333333333*(pow(x3,4)*(5*pow(x2,2) - 4*x2*x3 + pow(x3,2))*z)/(pow(x1,2)*pow(x1 - x2,2)*pow(x1 - x3,2));
				matVne[6] = (pow(x3,4)*(5*x1*x2 - 2*x1*x3 - 2*x2*x3 + pow(x3,2))*z)/(30.*x1*pow(x1 - x2,2)*x2*(x1 - x3)*(x2 - x3));
				matVne[7] = matVne[7] = (pow(x3,2)*(x1*(-10*pow(x2,2) + 10*x2*x3 - 3*pow(x3,2)) + x3*(5*pow(x2,2) - 6*x2*x3 + 2*pow(x3,2)))*z)/(30.*x1*(x1 - x2)*pow(x1 - x3,2)*(x2 - x3));
				matVne[8] = matVne[2];
				matVne[9] = (pow(x3,4)*(5*x1*x2 - 2*x1*x3 - 2*x2*x3 + pow(x3,2))*z)/(30.*x1*pow(x1 - x2,2)*x2*(x1 - x3)*(x2 - x3));
				matVne[10] = -0.03333333333333333*(pow(x3,4)*(5*pow(x1,2) - 4*x1*x3 + pow(x3,2))*z)/(pow(x1 - x2,2)*pow(x2,2)*pow(x2 - x3,2));
				matVne[11] = (pow(x3,2)*(2*x1*(5*x2 - 3*x3)*x3 + 5*pow(x1,2)*(-2*x2 + x3) + pow(x3,2)*(-3*x2 + 2*x3))*z)/(30.*x2*(-x1 + x2)*(x1 - x3)*pow(x2 - x3,2));
				matVne[12] = 0.0;
				matVne[13] = matVne[7];
				matVne[14] = matVne[11];
				matVne[15] = ((pow(x3,2)*(-15*pow(x2,2) + 24*x2*x3 - 10*pow(x3,2)) - 5*pow(x1,2)*(6*pow(x2,2) - 8*x2*x3 + 3*pow(x3,2)) + 4*x1*x3*(10*pow(x2,2) - 15*x2*x3 + 6*pow(x3,2)))*z)/(30.*pow(x1 - x3,2)*pow(x2 - x3,2));


			}
			else{
				matVne[0] = (z*(10*pow(x0,6) - 24*pow(x0,5)*(x1 + x2 + x3) + pow(x3,2)*(90*pow(x1,2)*pow(x2,2) + 40*x1*x2*(x1 + x2)*x3 - 5*(pow(x1,2) + 4*x1*x2 + pow(x2,2))*pow(x3,2) + 4*(x1 + x2)*pow(x3,3) - pow(x3,4)) - 120*x0*x1*x2*x3*(x2*x3 + x1*(x2 + x3)) + 15*pow(x0,4)*(pow(x1,2) + pow(x2,2) + 4*x2*x3 + pow(x3,2) + 4*x1*(x2 + x3)) - 40*pow(x0,3)*(pow(x1,2)*(x2 + x3) + x2*x3*(x2 + x3) + x1*(pow(x2,2) + 4*x2*x3 + pow(x3,2))) + 30*pow(x0,2)*(pow(x2,2)*pow(x3,2) + 4*x1*x2*x3*(x2 + x3) + pow(x1,2)*(pow(x2,2) + 4*x2*x3 + pow(x3,2))) + 60*pow(x1,2)*pow(x2,2)*pow(x3,2)*(log(x0) - log(x3))))/(30.*pow(x0 - x1,2)*pow(x0 - x2,2)*pow(x0 - x3,2));
				matVne[1] = (z*(2*pow(x0,6) + x0*pow(x3,2)*(-30*x1*pow(x2,2) - 20*x2*(2*x1 + x2)*x3 + 5*(x1 + 2*x2)*pow(x3,2) - 2*pow(x3,3)) + pow(x3,3)*(-20*x1*pow(x2,2) + 5*x2*(2*x1 + x2)*x3 - 2*(x1 + 2*x2)*pow(x3,2) + pow(x3,3)) - 3*pow(x0,5)*(x1 + 2*(x2 + x3)) + 30*pow(x0,2)*x2*x3*(x2*x3 + 2*x1*(x2 + x3)) + 5*pow(x0,4)*(pow(x2,2) + 4*x2*x3 + pow(x3,2) + 2*x1*(x2 + x3)) - 10*pow(x0,3)*(2*x2*x3*(x2 + x3) + x1*(pow(x2,2) + 4*x2*x3 + pow(x3,2))) + 60*x0*x1*pow(x2,2)*pow(x3,2)*(-log(x0) + log(x3))))/(30.*pow(x0 - x1,2)*(x0 - x2)*(x1 - x2)*(x0 - x3)*(x1 - x3));
				matVne[2] = (z*(-2*pow(x0,6) + 3*pow(x0,5)*(2*x1 + x2 + 2*x3) + pow(x3,3)*(20*pow(x1,2)*x2 - 5*x1*(x1 + 2*x2)*x3 + 2*(2*x1 + x2)*pow(x3,2) - pow(x3,3)) + x0*pow(x3,2)*(30*pow(x1,2)*x2 + 20*x1*(x1 + 2*x2)*x3 - 5*(2*x1 + x2)*pow(x3,2) + 2*pow(x3,3)) - 30*pow(x0,2)*x1*x3*(2*x2*x3 + x1*(2*x2 + x3)) - 5*pow(x0,4)*(pow(x1,2) + x3*(2*x2 + x3) + 2*x1*(x2 + 2*x3)) + 10*pow(x0,3)*(x2*pow(x3,2) + 2*x1*x3*(2*x2 + x3) + pow(x1,2)*(x2 + 2*x3)) + 60*x0*pow(x1,2)*x2*pow(x3,2)*(log(x0) - log(x3))))/(30.*(x0 - x1)*pow(x0 - x2,2)*(x1 - x2)*(x0 - x3)*(x2 - x3));
				matVne[3] = (z*(pow(x0,2)*(2*pow(x0,4) + 30*pow(x1,2)*pow(x2,2) - 6*pow(x0,3)*(x1 + x2) - 20*x0*x1*x2*(x1 + x2) + 5*pow(x0,2)*(pow(x1,2) + 4*x1*x2 + pow(x2,2))) + pow(x0,2)*(-3*pow(x0,3) + 10*pow(x0,2)*(x1 + x2) + 60*x1*x2*(x1 + x2) - 10*x0*(pow(x1,2) + 4*x1*x2 + pow(x2,2)))*x3 - 30*x1*x2*(x1*x2 + 2*x0*(x1 + x2))*pow(x3,2) + 10*(2*x1*x2*(x1 + x2) + x0*(pow(x1,2) + 4*x1*x2 + pow(x2,2)))*pow(x3,3) - 5*(pow(x1,2) + 4*x1*x2 + pow(x2,2) + 2*x0*(x1 + x2))*pow(x3,4) + 3*(x0 + 2*(x1 + x2))*pow(x3,5) - 2*pow(x3,6) + 60*x0*pow(x1,2)*pow(x2,2)*x3*(-log(x0) + log(x3))))/(30.*(x0 - x1)*(x0 - x2)*pow(x0 - x3,2)*(x1 - x3)*(x2 - x3));
				matVne[4] = matVne[1];
				matVne[5] = (z*(pow(x0,4)*(pow(x0,2) - 4*x0*x2 + 5*pow(x2,2)) - 4*pow(x0,3)*(pow(x0,2) - 5*x0*x2 + 10*pow(x2,2))*x3 + 5*pow(x0,3)*(x0 - 8*x2)*pow(x3,2) + 40*x0*x2*(x0 + x2)*pow(x3,3) - 5*(pow(x0,2) + 4*x0*x2 + pow(x2,2))*pow(x3,4) + 4*(x0 + x2)*pow(x3,5) - pow(x3,6) + 60*pow(x0,2)*pow(x2,2)*pow(x3,2)*(log(x0) - log(x3))))/(30.*pow(x0 - x1,2)*pow(x1 - x2,2)*pow(x1 - x3,2));
				matVne[6] = (z*(-pow(x0,6) + 5*pow(x0,2)*pow(x3,3)*(-4*(x1 + x2) + x3) + 2*pow(x0,5)*(x1 + x2 + 2*x3) + 20*pow(x0,3)*x3*(2*x1*x2 + (x1 + x2)*x3) + 2*x0*pow(x3,3)*(-20*x1*x2 + 5*(x1 + x2)*x3 - 2*pow(x3,2)) + pow(x3,4)*(5*x1*x2 - 2*(x1 + x2)*x3 + pow(x3,2)) - 5*pow(x0,4)*(x1*x2 + 2*(x1 + x2)*x3 + pow(x3,2)) + 60*pow(x0,2)*x1*x2*pow(x3,2)*(-log(x0) + log(x3))))/(30.*(x0 - x1)*(x0 - x2)*pow(x1 - x2,2)*(x1 - x3)*(x2 - x3));
				matVne[7] = (z*(pow(x0,6) - 2*pow(x0,5)*(x1 + 2*x2 + x3) + 2*x0*pow(x3,2)*(30*x1*pow(x2,2) - 10*x2*(2*x1 + x2)*x3 + 5*(x1 + 2*x2)*pow(x3,2) - 3*pow(x3,3)) + 5*pow(x0,2)*x3*(-6*x1*pow(x2,2) + 6*x2*(2*x1 + x2)*x3 - 2*(x1 + 2*x2)*pow(x3,2) + pow(x3,3)) + pow(x3,3)*(-10*x1*pow(x2,2) + 5*x2*(2*x1 + x2)*x3 - 3*(x1 + 2*x2)*pow(x3,2) + 2*pow(x3,3)) - 20*pow(x0,3)*x2*(x2*x3 + x1*(x2 + 2*x3)) + 5*pow(x0,4)*(x1*(2*x2 + x3) + x2*(x2 + 2*x3)) + 60*pow(x0,2)*x1*pow(x2,2)*x3*(log(x0) - log(x3))))/(30.*(x0 - x1)*(x1 - x2)*(x0 - x3)*pow(x1 - x3,2)*(x2 - x3));
				matVne[8] = matVne[2];
				matVne[9] = matVne[6];
				matVne[10] = (z*(pow(x0,4)*(pow(x0,2) - 4*x0*x1 + 5*pow(x1,2)) - 4*pow(x0,3)*(pow(x0,2) - 5*x0*x1 + 10*pow(x1,2))*x3 + 5*pow(x0,3)*(x0 - 8*x1)*pow(x3,2) + 40*x0*x1*(x0 + x1)*pow(x3,3) - 5*(pow(x0,2) + 4*x0*x1 + pow(x1,2))*pow(x3,4) + 4*(x0 + x1)*pow(x3,5) - pow(x3,6) + 60*pow(x0,2)*pow(x1,2)*pow(x3,2)*(log(x0) - log(x3))))/(30.*pow(x0 - x2,2)*pow(x1 - x2,2)*pow(x2 - x3,2));
				matVne[11] = (z*(pow(x0,3)*(pow(x0,3) - 20*pow(x1,2)*x2 - 2*pow(x0,2)*(2*x1 + x2) + 5*x0*x1*(x1 + 2*x2)) - pow(x0,2)*(2*pow(x0,3) + 30*pow(x1,2)*x2 - 5*pow(x0,2)*(2*x1 + x2) + 20*x0*x1*(x1 + 2*x2))*x3 + 30*x0*x1*(x0*x1 + 2*(x0 + x1)*x2)*pow(x3,2) - 10*(2*x0*x1*(x0 + x1) + (pow(x0,2) + 4*x0*x1 + pow(x1,2))*x2)*pow(x3,3) + 5*(pow(x0,2) + 4*x0*x1 + pow(x1,2) + 2*(x0 + x1)*x2)*pow(x3,4) - 3*(2*(x0 + x1) + x2)*pow(x3,5) + 2*pow(x3,6) + 60*pow(x0,2)*pow(x1,2)*x2*x3*(log(x0) - log(x3))))/(30.*(x0 - x2)*(-x1 + x2)*(x0 - x3)*(x1 - x3)*pow(x2 - x3,2));
				matVne[12] = matVne[3];
				matVne[13] = matVne[7];
				matVne[14] = matVne[11];
				matVne[15] = matVne[15] = (z*(pow(x0,2)*(pow(x0,4) - 90*pow(x1,2)*pow(x2,2) - 4*pow(x0,3)*(x1 + x2) - 40*x0*x1*x2*(x1 + x2) + 5*pow(x0,2)*(pow(x1,2) + 4*x1*x2 + pow(x2,2))) + 120*x0*x1*x2*(x1*x2 + x0*(x1 + x2))*x3 - 30*(pow(x1,2)*pow(x2,2) + 4*x0*x1*x2*(x1 + x2) + pow(x0,2)*(pow(x1,2) + 4*x1*x2 + pow(x2,2)))*pow(x3,2) + 40*(pow(x0,2)*(x1 + x2) + x1*x2*(x1 + x2) + x0*(pow(x1,2) + 4*x1*x2 + pow(x2,2)))*pow(x3,3) - 15*(pow(x0,2) + pow(x1,2) + 4*x1*x2 + pow(x2,2) + 4*x0*(x1 + x2))*pow(x3,4) + 24*(x0 + x1 + x2)*pow(x3,5) - 10*pow(x3,6) + 60*pow(x0,2)*pow(x1,2)*pow(x2,2)*(log(x0) - log(x3))))/(30.*pow(x0 - x3,2)*pow(x1 - x3,2)*pow(x2 - x3,2));

			}
		}
	}


	return matVne;

}







#endif // _EXACT_Vr_INTEGRALS_H_




