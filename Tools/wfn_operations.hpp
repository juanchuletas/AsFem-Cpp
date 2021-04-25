#if !defined(_WFN_OPERATIONS)
#define _WFN_OPERATIONS
#include<iostream>
#include<cmath>
 void getWfnPhase(int nodes, int orb,int *phase,double *wfn);


#endif // _WFN_OPERATIONS





namespace asfem_tools{

    namespace wfn{

        void getWfnPhase(int nodes, int orb,int *phase,double *wfn){
			int i, answer; 
			double prev, val, der;
			prev = wfn[orb];
			i=1;
			answer = 1;
			while(i<nodes && answer!=0){
				val = wfn[i + orb*nodes];
				der = val-prev;
				if(fabs(der)>1E-5){
					answer=0;
				}
				prev = val;
				i++;
			}
			if(der>0){
				*phase = 1;
			}
			else{
				*phase = -1;
			}

		}
    }

}