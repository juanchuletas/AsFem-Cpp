#if !defined(_GRID_TOOLS_H_)
#define _GRID_TOOLS_H_

#include<iostream>
#include<cmath>
namespace grid_tools{
    class froese_fischer{
    
        public:
        static double rmf; //common value = 5;
        static double hmf; //common value = 1/32

        static double kernel(int i, double z);
        static double inverseKernel(double rN,int atomicN);
    };
}

#endif // _GRID_TOOLS_H_
