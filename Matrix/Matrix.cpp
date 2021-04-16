#include<iostream>
#include <iomanip>
#include "Matrix.hpp"

template<class T>
Matrix<T>::Matrix()
:matSize{0}{
    //std::cout<<"Default Constructor:\n";

    items = nullptr;
}
template<class T>
Matrix<T>::Matrix(int inputSize)
:matSize{inputSize}{
    //std::cout<<"Default Constructor: allocating memory\n";

    items = new T[matSize];

}
template<class T>
Matrix<T>::Matrix(int inputSize,T& data)
:matSize{inputSize}{
    //std::cout<<"Default Constructor: FILL WITH DATA MAT\n";

    items = new T[matSize];
    for(int i=0; i<matSize; i++){
        items[i] = data;
    }

}
template<class T>
Matrix<T>::Matrix(const Matrix<T>&source)
:matSize{source.matSize}{
    //std::cout<<"CPY CONSTRUCTOR\n";

    items = new T[matSize];
    for(int i=0; i<matSize; i++){
        items[i] = source.items[i];
    }

}
template<class T>
Matrix<T>::~Matrix(){
    //std::cout<<"MAT Destructor\n";

    delete []items;
}
//****************************
// *** METHODS *************
template<class T>
void Matrix<T>::fillWithZeros(){
    for(int i=0; i<matSize; i++){
        items[i] = static_cast<T>(0);
    }
}
template<class T>
void Matrix<T>::printMatrix(){
    for(int i=0; i<matSize; i++){
        std::cout<<std::setprecision(6)<<std::fixed<<items[i]<<std::endl;
    }
}
template<class T>
void Matrix<T>::printMatrix(int m, int n, std::string name){

    std::cout<<name<<std::endl;
    for( int i = 0; i < m; i++ ) {
                for( int j = 0; j < n; j++ ) printf( "%lf\t", items[i*n+j] );
                printf( "\n" );
        }

}
template<class T>
void Matrix<T>::sort(){
    int m=-1;
    static const int M=7, NSTACK=64;
	int i,ir,j,k,jstack=-1,l=0,n=matSize;
	double a;
	Matrix<int> istack(NSTACK);
	if (m>0) n = MIN(m,n);
	ir=n-1;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=items[j];
				for (i=j-1;i>=l;i--) {
					if (items[i] <= a) break;
					items[i+1]=items[i];
				}
				items[i+1]=a;
			}
			if (jstack < 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(items[k],items[l+1]);
			if (items[l] > items[ir]) {
				SWAP(items[l],items[ir]);
			}
			if (items[l+1] > items[ir]) {
				SWAP(items[l+1],items[ir]);
			}
			if (items[l] > items[l+1]) {
				SWAP(items[l],items[l+1]);
			}
			i=l+1;
			j=ir;
			a=items[l+1];
			for (;;) {
				do i++; while (items[i] < a);
				do j--; while (items[j] > a);
				if (j < i) break;
				SWAP(items[i],items[j]);
			}
			items[l+1]=items[j];
			items[j]=a;
			jstack += 2;
			if (jstack >= NSTACK) throw("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}
//*** END METHODS **********

//******* GETTERS ***********
template<class T>
int Matrix<T>::getSize() const{
    return matSize;
}
//*** END GETTERS ***********
// ** SETTERS
template<class T>
void Matrix<T>::setMatrix(int inputSize){
    matSize = inputSize;
    items = new T[matSize];

}
// ** OPERATORS **
template<class T>
const T& Matrix<T>::operator[](int i)const{
    try{
        if (i<0 || i>=matSize) {
            throw 0;
        }
	    return items[i];
    }
    catch(int &ex){

        std::cerr<<"OUT OF BOUNDS\n"<<std::endl;
    }
}
template<class T>
T& Matrix<T>::operator[](int i){
    try{
        if (i<0 || i>=matSize) {
            throw 0;
        }
	    return items[i];
    }
    catch(int &ex){

        std::cerr<<"OUT OF BOUNDS\n"<<std::endl;
    }
}
template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &source){
    if(this!=&source){
        if(this->matSize==source.matSize){
            for(int i=0; i<matSize; i++){
                items[i] = source.items[i];
            }
        }
        else{
            std::cout<<"MATRICES HAVE DIFFERENT SIZES::FILLING WITH ZEROES\n";
            this->fillWithZeros();
        }
    }

}
// **END OPERATORS**