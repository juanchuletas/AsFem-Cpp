#if !defined(_MATRIX_H_)
#define _MATRIX_H_
#include<string>
template<class T> class Matrix{

    int matSize;
    T  *items;

    public:
        Matrix();
        Matrix(int inputSize);//Default constructor
        Matrix(int inputSize,T& data);//Copy Constructor
        Matrix(const Matrix &source);
        ~Matrix();
        //METHODS
        void fillWithZeros();
        void printMatrix();
        void printMatrix(int, int, std::string name);
        //GETTERS
        int getSize() const;
        //SETTERS
        void setMatrix(int dimension);

        //OPERATORS
        Matrix<T>& operator=(const Matrix<T>&source);
        const T& operator[](int i) const;
        T& operator[](int i);
};

#endif // _MATRIX_H_
