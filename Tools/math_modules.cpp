#include "math_modules.hpp"
void poissonSolver(double *hartree_vec,double *uij,double *lij,double *rho,int pnodes,double *b,double hp)
{
	
	int r_nodes = pnodes;
	int totSize = r_nodes*r_nodes;
	int N=r_nodes;
	int NRHS=1;
	int LDA=N;
	int LDB=N;
	int ipiv[N];
	int info;
	int Norb=1;
	double sum=0.0;
	double coeffwfn,coeffrho;
	const double PI = 3.14159265358979323846;
  //double *right_vec = new double[r_nodes];
	//LOCAL ARRAYS:
    double *right_vec = new double[r_nodes];
    double *aux_vec = new double[r_nodes];
    double *aux_mat = new double[r_nodes*r_nodes];
		
  ScalarXMatrix(1.0,lij,aux_mat,r_nodes,r_nodes);
		
  
	MatrixProduct(uij,rho,right_vec,r_nodes,r_nodes,NRHS);
  
	 //MatXvec(sij,rho,right_vec,r_nodes);
	//double qtot = 0.0;
	for(int i=0; i<r_nodes; i++){
		aux_vec[i] = right_vec[i] - b[i]*hp;
		//qtot  = qtot + right_vec[i];
		//printf("bvec[%d] = %lf\n",i,b[i]);
	}
	//printf("THE SUM IS: %lf\n",qtot);


	dgesv_(&N,&NRHS,aux_mat,&LDA,ipiv,aux_vec,&LDB,&info);
	
	

	if( info > 0 ) 
	{
        	printf( "The diagonal element of the triangular factor of A,\n" );
       		printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
        	printf( "the solution could not be computed.\n" );
        	exit(1);
	}

	for(int i=0; i<r_nodes; i++)
	{
		hartree_vec[i] = aux_vec[i];
	  //printf("HP[%d] = % lf\n",i,aux_vec[i]);
	}

	  delete [] aux_mat;
    delete [] right_vec;
    delete [] aux_vec;

}
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
