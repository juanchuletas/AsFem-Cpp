#if !defined(_GRID_H_)
#define _GRID_H_
#include<string>
#include <cmath>
#include "../Matrix/Matrix.hpp"
template<class T> class Grid {

    Matrix<T> grid;
    Matrix<T> elSize;
    int Ne;
    int order;
    int totalNodes;
    int polynomial;
    int bcNodes;
    int atomicN; //Froese-Fischer Grid needs it
    double r0,rN;
    std::string meshType;
    private://FUNCTIONS
        double FroeseFischer(int i, double z, double r, double h){return (exp(-r + h*(double)i)/z);}
        
    public:
        Grid();
        Grid(double rinit, double rfinal,int Ne, int order,std::string meshType);
         Grid(double rinit, double rfinal,int Ne, int order,std::string meshType, int atomicN);
        ~Grid();

        //METHODS

        void createGrid();
        void setGridData(double rinit, double rfinal,int Ne, int order,std::string meshType);
        void setGridData(double rinit, double rfinal,int Ne, int order,std::string meshType,int atomicN);
        void buildAtomic(int amtomicN);
        void buildChebyshev();
        T& getElementSize(int i);
        T& getGridItem(int i);
        int size();
        void printGrid();
        T& operator[](int i);

        
};
#include "Grid.cpp"
#endif // _GRID_H_
