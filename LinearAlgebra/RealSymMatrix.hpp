#if !defined(_LINEAR_ALGEBRA_H_)
#define _LINEAR_ALGEBRA_H_
#include <cmath>
#include <limits>       // std::numeric_limits
#include "../Matrix/Matrix.hpp"
#include "inlineStuff.hpp"
class RealSymMat{

    int Ntot, Mtot;
    Matrix<double> z;
    Matrix<double> d,e; 
    bool yVecs;
/*     inline  double SIGN(const double &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
    inline double SQR(const double a) {return a*a;} */
    void eigsrt(bool vec);

    public:
        RealSymMat(const Matrix<double> &targetMat,int Nrows,int Mcols);
        ~RealSymMat();
        void getTrid();
        void getQL();
        double pythag(double a,double b);
        void printResult();
        //void sort(Matrix<double> &arr, int);


};
RealSymMat::RealSymMat(const Matrix<double> &targetMat,int Nrows,int Mcols)
:z{targetMat},Ntot{Nrows}{
    yVecs = true;
    e.setMatrix(Ntot);
    d.setMatrix(Ntot);
    getTrid();
    getQL();
    eigsrt(true);
}
void RealSymMat::getTrid(){
    double h,scale,doubH,g,f;

    for (int i=Ntot-1;i>0;i--) {
		int l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (int k=0;k<i;k++)
				scale += abs(z[i*Ntot+k]);
			if (scale == 0.0)
				e[i]=z[i*Ntot + l];
			else {
				for (int k=0;k<i;k++) {
					z[i*Ntot + k] /= scale;
					h += z[i*Ntot + k]*z[i*Ntot + k];
				}
				f=z[i*Ntot + l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				z[i*Ntot + l]=f-g;
				f=0.0;
				for (int j=0;j<i;j++) {
					if (yVecs)
						z[j*Ntot + i]=z[i*Ntot + j]/h;
					g=0.0;
					for (int k=0;k<j+1;k++)
						g += z[j*Ntot + k]*z[i*Ntot + k];
					for (int k=j+1;k<i;k++)
						g += z[k*Ntot + j]*z[i*Ntot + k];
					e[j]=g/h;
					f += e[j]*z[i*Ntot + j];
				}
				doubH=f/(h+h);
				for (int j=0;j<i;j++) {
					f=z[i*Ntot + j];
					e[j]=g=e[j]-doubH*f;
					for (int k=0;k<j+1;k++)
						z[j*Ntot + k] -= (f*e[k]+g*z[i*Ntot + k]);
				}
			}
		} else
			e[i]=z[i*Ntot + l];
		d[i]=h;
	}
	if (yVecs) d[0]=0.0;
	e[0]=0.0;
	for (int i=0;i<Ntot;i++) {
		if (yVecs) {
			if (d[i] != 0.0) {
				for (int j=0;j<i;j++) {
					g=0.0;
					for (int k=0;k<i;k++)
						g += z[i*Ntot + k]*z[k*Ntot + j];
					for (int k=0;k<i;k++)
						z[k*Ntot + j] -= g*z[k*Ntot + i];
				}
			}
			d[i]=z[i*Ntot + i];
			z[i*Ntot + i]=1.0;
			for (int j=0;j<i;j++) z[j*Ntot + i]=z[i*Ntot + j]=0.0;
		} else {
			d[i]=z[i*Ntot + i];
		}
	}


}
void RealSymMat::getQL(){
    
    //Int m,l,iter,i,k;
    int l, m,i, k;
	double s,r,p,g,f,dd,c,b;
	
    //const double EPS=std::numeric_limits<double>::epsilon();
    const double EPS = 0.0000000001;
    printf("EPS = %lf\n",EPS);
	for (int i=1;i<Ntot;i++) e[i-1]=e[i];
	e[Ntot-1]=0.0;
	for (l=0;l<Ntot;l++) {
		int iter=0;
		do {
			for (m=l;m<Ntot-1;m++) {
				dd=abs(d[m])+abs(d[m+1]);
				if (abs(e[m]) <= EPS*dd) break;
			}
			if (m != l) {
				if (iter++ == 30) throw("Too many iterations in tqli");
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					if (yVecs) {
						for (k=0;k<Ntot;k++) {
							f=z[k*Ntot + (i+1)];
							z[k*Ntot + (i+1)]=s*z[k*Ntot + i]+c*f;
							z[k*Ntot + i]=c*z[k*Ntot + i]-s*f;
						}
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}

}
double RealSymMat::pythag(double a,double b){
    
    double absa=abs(a), absb=abs(b);
	return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
		(absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));

}

void RealSymMat::eigsrt(bool vec)
{
	int k;
	int n=d.getSize();
	for (int i=0;i<n-1;i++) {
		double p=d[k=i];
		for (int j=i;j<n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			if (vec ==true)
				for (int j=0;j<n;j++) {
					p=z[j*Ntot + i];
					z[j*Ntot + i]=z[j*Ntot + k];
					z[j*Ntot + k]=p;
				}
		}
	}
}
void RealSymMat::printResult(){


    std::cout<<"Results: \n";
    d.printMatrix(Ntot,Ntot,"Eigenvalues");
    z.printMatrix(Ntot,Ntot,"Eigenvectors");

}



// DESTRUCTOR
RealSymMat::~RealSymMat(){
    std::cout<<"Destructor\n";
}
#endif // _LINEAR_ALGEBRA_H_
 