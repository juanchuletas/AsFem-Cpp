#include"fixed_points_kinect_matrices.hpp"
double *getFixedPointsKinectMatrices(double *x,int order) {
	        
	int P = order +1;
    double *matK = new double[P*P];

	switch(order)
	{
		case 2:
		{
			matK[0] = 1/(-x[0] + x[2]) + (-x[0] + x[2])/(3.*pow(x[0] - x[1],2));
        	matK[1] = pow(x[0] - x[2],2)/(3.*pow(x[0] - x[1],2)*(x[1] - x[2]));
        	matK[2] = (1/(-x[0] + x[1]) + 3/(x[0] - x[2]) + 1/(-x[1] + x[2]))/3.;
        	matK[3] = matK[1];
        	matK[4] = -0.3333333333333333*pow(x[0] - x[2],3)/(pow(x[0] - x[1],2)*pow(x[1] - x[2],2));
        	matK[5] = pow(x[0] - x[2],2)/(3.*(x[0] - x[1])*pow(x[1] - x[2],2));
        	matK[6] = matK[2];
        	matK[7] = matK[5];
        	matK[8] = 1/(-x[0] + x[2]) + (-x[0] + x[2])/(3.*pow(x[1] - x[2],2));
        	break;
		}
		case 3:
		{
			double x0,x1,x2,x3;
            x0 = x[0]; x1=x[1]; x2=x[2]; x3=x[3];
			matK[0] = ((-9*pow(x0,5))/5. + 3*pow(x0,4)*(x1 + x2 + x3) - x0*pow(x2*x3 + x1*(x2 + x3),2) - (2*pow(x0,3)*(2*pow(x1,2) + 2*pow(x2,2) + 7*x2*x3 + 2*pow(x3,2) + 7*x1*(x2 + x3)))/3. + 2*pow(x0,2)*(pow(x1,2)*(x2 + x3) + x2*x3*(x2 + x3) + x1*(pow(x2,2) + 3*x2*x3 + pow(x3,2))) + (5*x1*(2*x2 - x3)*pow(x3,3) + pow(x3,3)*(5*pow(x2,2) - 5*x2*x3 + 2*pow(x3,2)) + 5*pow(x1,2)*(3*pow(x2,2)*x3 + pow(x3,3)))/15.)/(pow(x0 - x1,2)*pow(x0 - x2,2)*pow(x0 - x3,2));
			matK[1] = (pow(x0 - x3,2)*(9*pow(x0,2) + 10*x2*(x1 + x2) - 5*x0*(x1 + 4*x2) + 7*x0*x3 - 5*(x1 + 2*x2)*x3 + 4*pow(x3,2)))/(30.*pow(x0 - x1,2)*(x0 - x2)*(x1 - x2)*(x1 - x3));
			matK[2] = -0.03333333333333333*(pow(x0 - x3,2)*(9*pow(x0,2) + 10*x1*(x1 + x2) - 5*x0*(4*x1 + x2) + 7*x0*x3 - 5*(2*x1 + x2)*x3 + 4*pow(x3,2)))/((x0 - x1)*pow(x0 - x2,2)*(x1 - x2)*(x2 - x3));
			matK[3] = (9*pow(x0,4) + 30*pow(x1,2)*pow(x2,2) - 30*x1*x2*(x1 + x2)*x3 + 10*(pow(x1,2) + 5*x1*x2 + pow(x2,2))*pow(x3,2) - 20*(x1 + x2)*pow(x3,3) + 9*pow(x3,4) + 4*pow(x0,3)*(-5*(x1 + x2) + x3) + 2*pow(x0,2)*(5*(pow(x1,2) + 5*x1*x2 + pow(x2,2)) - 5*(x1 + x2)*x3 + 2*pow(x3,2)) + 2*x0*(-15*x1*x2*(x1 + x2) + 5*pow(x1 + x2,2)*x3 - 5*(x1 + x2)*pow(x3,2) + 2*pow(x3,3)))/(30.*(x0 - x1)*(x0 - x2)*(x0 - x3)*(x1 - x3)*(x2 - x3));
			matK[4] = matK[1];
			matK[5] = -0.06666666666666667*(pow(x0 - x3,3)*(2*pow(x0,2) + 5*pow(x2,2) - 5*x2*x3 + 2*pow(x3,2) + x0*(-5*x2 + x3)))/(pow(x0 - x1,2)*pow(x1 - x2,2)*pow(x1 - x3,2));matK[6] = (pow(x0 - x3,3)*(4*pow(x0,2) + 10*x1*x2 - 5*x0*(x1 + x2) + 2*x0*x3 - 5*(x1 + x2)*x3 + 4*pow(x3,2)))/(30.*(x0 - x1)*(x0 - x2)*pow(x1 - x2,2)*(x1 - x3)*(x2 - x3));
			matK[6] = (pow(x0 - x3,3)*(4*pow(x0,2) + 10*x1*x2 - 5*x0*(x1 + x2) + 2*x0*x3 - 5*(x1 + x2)*x3 + 4*pow(x3,2)))/(30.*(x0 - x1)*(x0 - x2)*pow(x1 - x2,2)*(x1 - x3)*(x2 - x3));
			matK[7] = -0.03333333333333333*(pow(x0 - x3,2)*(4*pow(x0,2) + 10*x2*(x1 + x2) - 5*x0*(x1 + 2*x2) + 7*x0*x3 - 5*(x1 + 4*x2)*x3 + 9*pow(x3,2)))/((x0 - x1)*(x1 - x2)*pow(x1 - x3,2)*(x2 - x3));
			matK[8] = matK[2];
            matK[9] = matK[6];
			matK[10] = -0.06666666666666667*(pow(x0 - x3,3)*(2*pow(x0,2) + 5*pow(x1,2) - 5*x1*x3 + 2*pow(x3,2) + x0*(-5*x1 + x3)))/(pow(x0 - x2,2)*pow(x1 - x2,2)*pow(x2 - x3,2));
			matK[11] = -0.03333333333333333*(pow(x0 - x3,2)*(4*pow(x0,2) + 10*x1*(x1 + x2) - 5*x0*(2*x1 + x2) + 7*x0*x3 - 5*(4*x1 + x2)*x3 + 9*pow(x3,2)))/((x0 - x2)*(-x1 + x2)*(x1 - x3)*pow(x2 - x3,2));
			matK[12] = matK[3];
          	matK[13] = matK[7];
            matK[14] = matK[11];
			matK[15] = (-2*pow(x0,4) - 15*pow(x1,2)*pow(x2,2) + pow(x0,3)*(5*(x1 + x2) - 2*x3) + 30*x1*x2*(x1 + x2)*x3 - 10*(2*pow(x1,2) + 7*x1*x2 + 2*pow(x2,2))*pow(x3,2) + 45*(x1 + x2)*pow(x3,3) - 27*pow(x3,4) + pow(x0,2)*(-5*pow(x1 + x2,2) + 5*(x1 + x2)*x3 - 2*pow(x3,2)) + x0*x3*(10*pow(x1 + x2,2) - 25*(x1 + x2)*x3 + 18*pow(x3,2)))/(15.*(x0 - x3)*pow(x1 - x3,2)*pow(x2 - x3,2));
			break;

		}
	}

return matK;

}