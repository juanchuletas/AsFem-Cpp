//This is the main file to change stuff about the functions to build the 
//kernels
#include "grid_tools.hpp"

namespace grid_tools{

    double froese_fischer::kernel(int i, double z){
        double x,kern;
        kern = exp(-rmf+hmf*(double)i); 
        x = kern/z;
        return x;
    }
    double froese_fischer::inverseKernel(double rN,int atomicN){
        double z = static_cast<double>(atomicN);
        double x = (1.0/hmf)*(log(rN*z) + rmf);
        return x;
    }
}
double grid_tools::froese_fischer::hmf;
double grid_tools::froese_fischer::rmf;