#include "fixed_elemets_kinect_matrices.hpp"
double *ElementalKinectMatrix(int order)
{
        int P = order + 1;

        double *eMatK = new double[P*P];

        switch(order)
        {
                case 1:

                        eMatK[0] = 0.5; eMatK[1] = -0.5;
                        eMatK[2] = eMatK[1]; eMatK[3] = eMatK[0];
                        break;
                case 2:
                        eMatK[0] = 7.0/6.0;      eMatK[1] = -4.0/3.0;    eMatK[2] = 1.0/6.0;
                        eMatK[3] = eMatK[1];     eMatK[4] = 8.0/3.0;     eMatK[5] = -4.0/3.0;
                        eMatK[6] = eMatK[2];     eMatK[7] = eMatK[5];     eMatK[8] = 7.0/6.0;
                        break;
                case 3:
                        eMatK[0] = 37./20.;        eMatK[1] = -189./80.;   eMatK[2] = 27./40.;    eMatK[3] = -13./80.;
                        eMatK[4] = eMatK[1];       eMatK[5] = 27./5.;      eMatK[6] = -297./80.;  eMatK[7] = 27./40.;
                        eMatK[8] = eMatK[2];       eMatK[9] = eMatK[6];       eMatK[10]= 27./5.;     eMatK[11]= -189./80.;
                        eMatK[12]= eMatK[3];       eMatK[13]= eMatK[7];       eMatK[14]= eMatK[11];     eMatK[15]= 37./20.;
                        break;
                case 4:

                        eMatK[0]  =  2.6058201058201058201;     eMatK[1]  = -3.6232804232804232804;
                        eMatK[2]  =  1.6126984126984126984;     eMatK[3]  = -0.77883597883597883598;
                        eMatK[4]  =  0.18359788359788359788;    eMatK[5]  = -3.6232804232804232804;
                        eMatK[6]  =  8.8042328042328042328;     eMatK[7]  = -7.5174603174603174603;
                        eMatK[8]  =  3.1153439153439153439;     eMatK[9]  = -0.77883597883597883598;
                        eMatK[10] =  1.6126984126984126984;     eMatK[11] = -7.5174603174603174603;
                        eMatK[12] = 11.809523809523809524;      eMatK[13] = -7.5174603174603174603;
                        eMatK[14] =  1.6126984126984126984;     eMatK[15] = -0.77883597883597883598;
                        eMatK[16] =  3.1153439153439153439;     eMatK[17] = -7.5174603174603174603;
                        eMatK[18] =  8.8042328042328042328;     eMatK[19] = -3.6232804232804232804;
                        eMatK[20] =  0.18359788359788359788;    eMatK[21] = -0.77883597883597883598;
                        eMatK[22] =  1.6126984126984126984;     eMatK[23] = -3.6232804232804232804;
                        eMatK[24] =  2.6058201058201058201;
                        break;
        }

        return eMatK;
}
