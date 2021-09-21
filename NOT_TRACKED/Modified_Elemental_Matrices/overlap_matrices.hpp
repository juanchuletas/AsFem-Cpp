#if !defined(_MOD_OVERLAP_MATRICES_H_)
#define _MOD_OVERLAP_MATRICES_H_

#include<iostream>
#include<cmath>

double *getOverlapElementalMatrices(double *x,int order){
    int P = order + 1;

	double *matS = new double[P*P];
    switch(order)
    {
        case 2:
        {
            double x0,x1,x2;
            x0 = x[0]; x1 = x[1]; x2 = x[2];
            matS[0] = -0.03333333333333333*((x0 - x2)*(6*pow(x0,2) + 10*pow(x1,2) - 5*x1*x2 + pow(x2,2) + 3*x0*(-5*x1 + x2)))/pow(x0 - x1,2);
            matS[1]  = -0.016666666666666666*(pow(x0 - x2,3)*(3*x0 - 5*x1 + 2*x2))/(pow(x0 - x1,2)*(x1 - x2));
            matS[2] = ((x0 - x2)*(3*pow(x0,2) - 10*x0*x1 + 10*pow(x1,2) + 4*x0*x2 - 10*x1*x2 + 3*pow(x2,2)))/(60.*(x0 - x1)*(x1 - x2));
            matS[3] = matS[1];
            matS[4] = -0.03333333333333333*pow(x0 - x2,5)/(pow(x0 - x1,2)*pow(x1 - x2,2));
            matS[5] = (pow(x0 - x2,3)*(2*x0 - 5*x1 + 3*x2))/(60.*(x0 - x1)*pow(x1 - x2,2));
            matS[6] = matS[2];
            matS[7] = matS[5];
            matS[8] = -0.03333333333333333*((x0 - x2)*(pow(x0,2) - 5*x0*x1 + 10*pow(x1,2) + 3*x0*x2 - 15*x1*x2 + 6*pow(x2,2)))/pow(x1 - x2,2);
            break;

        }
        case 3:
        {
            double x0,x1,x2,x3;
            x0 = x[0]; x1=x[1]; x2=x[2]; x3=x[3];
            matS[0] = -0.004761904761904762*((x0 - x3)*(30*pow(x0,4) + 70*pow(x1,2)*pow(x2,2) - 10*pow(x0,3)*(7*(x1 + x2) - 2*x3) - 35*x1*x2*(x1 + x2)*x3 + 7*(pow(x1,2) + 4*x1*x2 + pow(x2,2))*pow(x3,2) - 7*(x1 + x2)*pow(x3,3) + 2*pow(x3,4) + 6*pow(x0,2)*(7*(pow(x1,2) + 4*x1*x2 + pow(x2,2)) - 7*(x1 + x2)*x3 + 2*pow(x3,2)) + 3*x0*(-35*x1*x2*(x1 + x2) + 7*(pow(x1,2) + 4*x1*x2 + pow(x2,2))*x3 - 7*(x1 + x2)*pow(x3,2) + 2*pow(x3,3))))/(pow(x0 - x1,2)*pow(x0 - x2,2));
            matS[1] = -0.002380952380952381*(pow(x0 - x3,3)*(10*pow(x0,3) - 35*x1*pow(x2,2) - 2*pow(x0,2)*(7*x1 + 14*x2 - 6*x3) + 14*x2*(2*x1 + x2)*x3 - 7*(x1 + 2*x2)*pow(x3,2) + 4*pow(x3,3) + x0*(21*x2*(2*x1 + x2) - 14*(x1 + 2*x2)*x3 + 9*pow(x3,2))))/(pow(x0 - x1,2)*(x0 - x2)*(x1 - x2)*(x1 - x3));
            matS[2] = (pow(x0 - x3,3)*(10*pow(x0,3) - 35*pow(x1,2)*x2 - 2*pow(x0,2)*(14*x1 + 7*x2 - 6*x3) + 14*x1*(x1 + 2*x2)*x3 - 7*(2*x1 + x2)*pow(x3,2) + 4*pow(x3,3) + x0*(21*x1*(x1 + 2*x2) - 14*(2*x1 + x2)*x3 + 9*pow(x3,2))))/(420.*(x0 - x1)*pow(x0 - x2,2)*(x1 - x2)*(x2 - x3));
            matS[3] = -0.002380952380952381*((x0 - x3)*(10*pow(x0,4) + 70*pow(x1,2)*pow(x2,2) - 4*pow(x0,3)*(7*(x1 + x2) - 4*x3) - 70*x1*x2*(x1 + x2)*x3 + 21*(pow(x1,2) + 4*x1*x2 + pow(x2,2))*pow(x3,2) - 28*(x1 + x2)*pow(x3,3) + 10*pow(x3,4) + 3*pow(x0,2)*(7*(pow(x1,2) + 4*x1*x2 + pow(x2,2)) - 14*(x1 + x2)*x3 + 6*pow(x3,2)) + 2*x0*(-35*x1*x2*(x1 + x2) + 14*(pow(x1,2) + 4*x1*x2 + pow(x2,2))*x3 - 21*(x1 + x2)*pow(x3,2) + 8*pow(x3,3))))/((x0 - x1)*(x0 - x2)*(x1 - x3)*(x2 - x3));
            matS[4] = matS[1];
            matS[5] = -0.004761904761904762*(pow(x0 - x3,5)*(2*pow(x0,2) - 7*x0*x2 + 7*pow(x2,2) + 3*x0*x3 - 7*x2*x3 + 2*pow(x3,2)))/(pow(x0 - x1,2)*pow(x1 - x2,2)*pow(x1 - x3,2));
            matS[6] = (pow(x0 - x3,5)*(4*pow(x0,2) + 14*x1*x2 - 7*x0*(x1 + x2) + 6*x0*x3 - 7*(x1 + x2)*x3 + 4*pow(x3,2)))/(420.*(x0 - x1)*(x0 - x2)*pow(x1 - x2,2)*(x1 - x3)*(x2 - x3));
            matS[7] = -0.002380952380952381*(pow(x0 - x3,3)*(4*pow(x0,3) - 35*x1*pow(x2,2) + 21*x2*(2*x1 + x2)*x3 - 14*(x1 + 2*x2)*pow(x3,2) + 10*pow(x3,3) + pow(x0,2)*(-7*(x1 + 2*x2) + 9*x3) + 2*x0*(7*x2*(2*x1 + x2) - 7*(x1 + 2*x2)*x3 + 6*pow(x3,2))))/((x0 - x1)*(x1 - x2)*pow(x1 - x3,2)*(x2 - x3));
            matS[8] = matS[2];
            matS[9] = matS[6];
            matS[10] = -0.004761904761904762*(pow(x0 - x3,5)*(2*pow(x0,2) - 7*x0*x1 + 7*pow(x1,2) + 3*x0*x3 - 7*x1*x3 + 2*pow(x3,2)))/(pow(x0 - x2,2)*pow(x1 - x2,2)*pow(x2 - x3,2));
            matS[11] = -0.002380952380952381*(pow(x0 - x3,3)*(4*pow(x0,3) - 35*pow(x1,2)*x2 + 21*x1*(x1 + 2*x2)*x3 - 14*(2*x1 + x2)*pow(x3,2) + 10*pow(x3,3) + pow(x0,2)*(-7*(2*x1 + x2) + 9*x3) + 2*x0*(7*x1*(x1 + 2*x2) - 7*(2*x1 + x2)*x3 + 6*pow(x3,2))))/((x0 - x2)*(-x1 + x2)*(x1 - x3)*pow(x2 - x3,2));
            matS[12] = matS[3];         
            matS[13] = matS[7];         
            matS[14] = matS[11];
            matS[15] = -0.004761904761904762*((x0 - x3)*(2*pow(x0,4) + 70*pow(x1,2)*pow(x2,2) - 105*x1*x2*(x1 + x2)*x3 + 42*(pow(x1,2) + 4*x1*x2 + pow(x2,2))*pow(x3,2) - 70*(x1 + x2)*pow(x3,3) + 30*pow(x3,4) + pow(x0,3)*(-7*(x1 + x2) + 6*x3) + pow(x0,2)*(7*(pow(x1,2) + 4*x1*x2 + pow(x2,2)) - 21*(x1 + x2)*x3 + 12*pow(x3,2)) + x0*(-35*x1*x2*(x1 + x2) + 21*(pow(x1,2) + 4*x1*x2 + pow(x2,2))*x3 - 42*(x1 + x2)*pow(x3,2) + 20*pow(x3,3))))/(pow(x1 - x3,2)*pow(x2 - x3,2));
            break;
      }

    }

return matS;

}



#endif // _MOD_OVERLAP_MATRICES_H_
