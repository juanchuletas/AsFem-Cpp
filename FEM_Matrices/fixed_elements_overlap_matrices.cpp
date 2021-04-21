#include"fixed_elements_overlap_matrices.hpp"
double *ElementalOverlapMatrix(int order)
{
        int P = order + 1;

        double *eMatS = new double[P*P];

        switch(order)
        {
                case 1:
                        eMatS[0] = 2.0/3.0;  eMatS[1] = 1.0/3.0;
                        eMatS[2] = eMatS[1]; eMatS[3] = eMatS[0];
                        break;
                case 2:

                        eMatS[0] = 4.0/15.0;       eMatS[1] = 2./15.0;       eMatS[2] = -1.0/15.0;
                        eMatS[3] = eMatS[1];     eMatS[4] = 16.0/15.0;      eMatS[5] = 2.0/15.0;
                        eMatS[6] = eMatS[2];     eMatS[7] = eMatS[5];     eMatS[8] = 4.0/15.0;
                        break;
                case 3:

                        eMatS[0] = 16./105.;       eMatS[1] = 33./280.;    eMatS[2] = -3./70.;      eMatS[3] = 19./840.;
                        eMatS[4] = eMatS[1];       eMatS[5] = 27./35.;     eMatS[6] = -27./280.;    eMatS[7] = -3./70.;
                        eMatS[8] = eMatS[2];       eMatS[9] = eMatS[6];    eMatS[10]= 27./35.;       eMatS[11]= 33./280.;
                        eMatS[12]= eMatS[3];       eMatS[13]= eMatS[7];    eMatS[14]= eMatS[11];     eMatS[15]= 16./105.;
                        break;
                case 4:
                        eMatS[0]  =  0.10299823633156966490;    eMatS[1]  =  0.10440917107583774250;
                        eMatS[2]  = -0.061375661375661375661;   eMatS[3]  =  0.019753086419753086420;
                        eMatS[4]  = -0.010229276895943562610;   eMatS[5]  =  0.10440917107583774250;
                        eMatS[6]  =  0.63209876543209876543;    eMatS[7]  = -0.13544973544973544974;
                        eMatS[8]  =  0.090299823633156966490;   eMatS[9]  =  0.019753086419753086420;
                        eMatS[10] = -0.061375661375661375661;   eMatS[11] = -0.13544973544973544974;
                        eMatS[12] =  0.66031746031746031746;    eMatS[13] = -0.13544973544973544974;
                        eMatS[14] = -0.061375661375661375661;   eMatS[15] =  0.019753086419753086420;
                        eMatS[16] =  0.090299823633156966490;   eMatS[17] = -0.13544973544973544974;
                        eMatS[18] =  0.63209876543209876543;    eMatS[19] =  0.10440917107583774250;
                        eMatS[20] = -0.010229276895943562610;   eMatS[21] =  0.019753086419753086420;
                        eMatS[22] = -0.061375661375661375661;   eMatS[23] =  0.10440917107583774250;
                        break;
        }

        return eMatS;
}
