#include "sycl_mathtools.hpp"


namespace SYCL_ENABLE{

    double * MathTools::rayleighQuotient(double *f_mat, double *smat, double *trialVec,double *exchangevec,int totOrbitals, int bcSize){
        
        printf("\x1B[32m***** Rayleigh Quotient SYCL Module *****\033[0m\t\t");
        printf("\n");
        //Creates  a Queue to link a gpu device
        auto Q = cl::sycl::queue{cl::sycl::gpu_selector{}};
         printf("\x1B[33mChosen device:\033[0m");
        std::cout << " "  << Q.get_device().get_info<cl::sycl::info::device::name>()<<std::endl;
        std::cout<< "Max Work Group Size: "<< Q.get_device().get_info<cl::sycl::info::device::max_work_group_size>()<<std::endl;

        double *matprodNum= new double[totOrbitals*bcSize];
        double *matprodDenom= new double[totOrbitals*bcSize];
        double *eigenValues = new double[totOrbitals];

        syclMatMul(Q, trialVec,totOrbitals,bcSize,f_mat,bcSize,bcSize,matprodNum);
        syclMatMul(Q,trialVec,totOrbitals,bcSize,smat,bcSize,bcSize,matprodDenom);
        double numerator,denominator;
        for(int i=0; i<(totOrbitals); i++){
            numerator = 0.0;
            denominator = 0.0;
            for(int j=0; j<bcSize; j++){
            
                numerator = numerator + trialVec[j + i*bcSize]*matprodNum[j + i*bcSize];
                denominator = denominator + trialVec[j + i*bcSize]*matprodDenom[j + i*bcSize];
            
            }
            eigenValues[i] = numerator/denominator;
            //printf("Orbital %d energy iterative SCF = %lf\n",i,0.5*orbitalEnergy[i]);
        } 
        



        delete [] matprodDenom;
        delete [] matprodNum;



        return eigenValues;
    }


}