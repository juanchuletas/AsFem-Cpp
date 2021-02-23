#include<iostream>
#include <iomanip>
#include "Matrix.hpp"

template<class T>
Matrix<T>::Matrix()
:matSize{0}{
    std::cout<<"Default Constructor:\n";

    items = nullptr;
}
template<class T>
Matrix<T>::Matrix(int inputSize)
:matSize{inputSize}{
    std::cout<<"Default Constructor: allocating memory\n";

    items = new T[matSize];

}
template<class T>
Matrix<T>::Matrix(int inputSize,T& data)
:matSize{inputSize}{
    std::cout<<"Default Constructor: FILL WITH DATA MAT\n";

    items = new T[matSize];
    for(int i=0; i<matSize; i++){
        items[i] = data;
    }

}
template<class T>
Matrix<T>::Matrix(const Matrix<T>&source)
:matSize{source.matSize}{
    std::cout<<"CPY CONSTRUCTOR\n";

    items = new T[matSize];
    for(int i=0; i<matSize; i++){
        items[i] = source.items[i];
    }

}
template<class T>
Matrix<T>::~Matrix(){
    std::cout<<"MAT Destructor\n";

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