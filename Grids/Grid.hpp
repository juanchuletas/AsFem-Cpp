#if !defined(_GRID_H_)
#define _GRID_H_
#include<string>
#include <cmath>
#include "../Matrix/Matrix.hpp"
#include "../Static/grid_tools.hpp"
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
        double FroeseFischer(int i, double z, double r, double h){
            double x = exp(-r + h*(double)i/z);
            return x;
        }
        
    public:
        Grid();
        Grid(double rinit, double rfinal,int Ne, int order,std::string meshType);
         Grid(double rinit, double rfinal,int Ne, int order,std::string meshType, int atomicN);
        ~Grid();

        //METHODS

        void createGrid();
        void createGrid(std::string name);
        void setSize(int _size);
        void setGridData(double rinit, double rfinal,int Ne, int order,std::string meshType);
        void setGridData(double rinit, double rfinal,int Ne, int order,std::string meshType,int atomicN);
        void exponential();
        void chebyshev();
        void froeseFischer(int atomicN);
        void buildAtomic(int amtomicN);
        void buildAtomic(int atomicN, std::string name);
        void buildChebyshev();
        int forcedInsertion(double cutVal);
        void refineNear(int points, double delta, int steps);
        T& getElementSize(int i);
        T& getGridItem(int i);
        int size();
        void printGrid();
        T& operator[](int i);
        Grid<T>& operator=(const Grid<T>&source);

        
};
#include "Grid.cpp"
#endif // _GRID_H_
