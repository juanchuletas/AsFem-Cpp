#include<iostream>
#include<cmath>
extern "C" void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
namespace asfem_tools{
	
	void diag (int n, double *h, double *s, double *e, double *v)
	{
	  int lwork, i, j, k, nn, ij, info;
	  static double small = 1.0e-10;
	  char name1[] = "V";
	  char name2[] = "U";
	  double *work, *aux, *uno;

	  lwork=3*n;
	  work = new double[lwork];
	  //work = (double *) malloc( lwork * sizeof (double));
		aux = new double[n*n];
	  //aux = (double *) malloc( n * n * sizeof (double));
	  for ( ij = 0; ij < n*n; ij++ )
	  {
		aux[ij] = s[ij];
	  }

	  dsyev_( name1, name2, &n, aux, &n, e, work, &lwork, &info );

	  if ( info !=0 )
	  {
		fprintf (stderr, "Overlap matrix diagonalization failed\n");
		exit (1);
	  }
	// print_matrix( "S matrix ", nbasis, nbasis, aux, nbasis );
	 //printf("\n");

	 nn = 0;
	  for ( i = 0; i < n; i++ )
	  {
		if ( e[i] > small)
		{
		  for ( j = 0; j < n; j++ )
		  {
			aux[j+nn*n] = aux[j+i*n] / sqrt(e[i]);
		  }
		  nn++;
		}
	  }
	  if ( nn < n )
	  {
		 fprintf (stdout, " # of linearly independent vectors = %d\n", nn);
	  }

	 //print_matrix( "X matrix ", nbasis, nbasis, aux, nbasis );
	 //printf("\n");

	  for ( i = 0; i < n; i++)
	  {
		for ( j = 0; j < nn; j++)
		{
		  v[i+j*n] = 0.0;
		  for  ( k = 0; k < n; k++)
		  {
			v[i+j*n] = v[i+j*n] + h[i+k*n] * aux[k+j*n];
		  }
		}
	  }
		uno = new double[nn*nn];
	  //uno = (double *) malloc( nn * nn * sizeof (double));
	  for ( i = 0; i < nn; i++ )
	  {
		for ( j = 0; j < nn; j++ )
		{
		  uno[i+j*nn] = 0.0;
		  for  ( k = 0; k < n; k++)
		  {
			uno[i+j*nn] = uno[i+j*nn] + aux[k+i*n] * v[k+j*n];
		  }
		}
	  }
	  //print_matrix( "F\' matrix ", nbasis, nbasis, uno, nbasis );
	  //printf("\n");
	   /* ---------------- Diagonalization of transformed F or H matrix -------------*/
	  dsyev_(name1,name2, &nn, uno, &nn, e, work, &lwork, &info );

	  if ( info !=0 )
	  {
		fprintf (stderr, "Fock matrix diagonalization failed\n");
		exit (1);
	  }


	  for ( i = 0; i < n; i++ )
	  {
		for ( j = 0; j < nn; j++ )
		{
		  v[i+j*n] = 0.0;
		  for  ( k = 0; k < n; k++ )
		  {
			v[i+j*n] = v[i+j*n] + aux[i+k*n] * uno[k+j*nn];
		  }
		}
	  }
		//print_matrix( "eigenvalues ", 1, nbasis, e, 1 );
		//printf("\n");
		delete [] uno;
		delete [] aux;
		delete [] work;
	  //free(uno); free(aux); free(work);
	  return;
	}

	double nPolExtrapolation(double *x,double *y, int N,double target)
	{
		double *c_vec,*d_vec,dy,ho,hp,den,w;
		double result,deltaTot,delta;
		dy = 0.00000000000001;
		int ns=0;
		c_vec = (double *)malloc(sizeof(double)*N);
		d_vec = (double *)malloc(sizeof(double)*N);
		delta = fabs(target-x[0]);
		for(int i=0; i<N; i++)
		{
			if((deltaTot=fabs(target-x[i])) < delta)
			{
				ns = i;
				delta = deltaTot;
			}
			c_vec[i] = y[i];
			d_vec[i] = y[i];
		}
		result = y[ns];
		ns = ns - 1;
		for(int m=0; m<N-1; m++)
		{
			for(int i=0; i<N-m-1; i++)
			{
				ho = x[i]-target;
				hp = x[i+m+1]-target;
				w = c_vec[i+1]-d_vec[i];
				if((den=ho-hp)==0.0){printf("Error in Extrapolation\n");exit(1);}
				den=w/den;
				d_vec[i] = hp*den;
				c_vec[i] = ho*den;
			}
			if( 2*(ns+1) < N-m-1 ) 
			{
				dy = c_vec[ns+1];
			} 
			else 
			{
				dy = d_vec[ns];
				ns = ns-1;
			}
			result = result + (dy);
		}
		return result;

	}
	double doInterpolation(double *xin,double *yin,int grade){
		
		double value;
		double x[8],y[8];	
		x[0]=xin[1];       y[0]=yin[1];
		x[1]=xin[2];       y[1]=yin[2];
		x[2]=xin[3];       y[2]=yin[3];
		x[3]=xin[4];       y[3]=yin[4];
		x[4]=xin[5];       y[4]=yin[5];
		x[5]=xin[6];       y[5]=yin[6];
		x[6]=xin[7];       y[6]=yin[7];
		x[7]=xin[8];       y[7]=yin[8];
		value = nPolExtrapolation(x,y,grade,0.000000);


		return value;
	}
	void sumArrays(double *h,double *v,double *k,int M,int N)
	{
		for(int i=0; i<M*N; i++)
		{
			h[i] = v[i] + k[i];
		}
	}
	int getFixedElements(std::string kindOfGrid,double rN, int order,int atomicN){
		int nele;
		double totpoints = atomicN*24.0*(log(rN) + 5);  
		int points = floor(totpoints);
		if(points%2==0){
			points = points - 1;
		}
		nele = (points-1)/order;
		return nele;
	}
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
	double getHFenergy(double *hcore, double *fmat, double *densMat,int totSize)
	{
		double energy = 0.0;

		for(int i=0; i<totSize; i++)
		{
			energy = energy + 0.5*(hcore[i]*densMat[i]);
			energy = energy + 0.5*(fmat[i]*densMat[i]);
			printf("e = %lf\n", energy);

		}
		return 0.5*energy;


	}




}//END OF NAMESPACE