#if !defined(_SYCL_MT_H_)
#define _SYCL_MT_H_
#include "../SYCL/SyclMatMul.hpp"
namespace SYCL_ENABLE{

    class MathTools{

        public:
            static double *rayleighQuotient(double *f_mat, double *smat, double *wfn,double *exchangevec,int totOrbitals, int bcSize);

    };

}

#endif // _SYCL_MT_H_
