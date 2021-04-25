//This is the main file to change stuff about the functions to build the 
//kernels
#include "grid_tools.hpp"

namespace grid_tools{

    double froese_fischer::kernel(int i, double z){
        double x; 
        x = exp(-rmf + hmf*(double)i/z);
    }
    double froese_fischer::inverseKernel(double rN,int atomicN){
        double z = static_cast<double>(atomicN);
        double x = z*(1.0/hmf)*(log(rN) + rmf);

    }
}
double grid_tools::froese_fischer::hmf;
double grid_tools::froese_fischer::rmf;