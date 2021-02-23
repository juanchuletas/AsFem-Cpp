#include<iostream>
#include<string>
//#include "UserInputData.cpp"

extern int read_data(char *input_name[],std::string &potential,std::string &atom, int &angular,int &Ne, double &r0, double &rc,int &order, std::string &mesh,double &lambda)
;
int main (int argc, char **argv){

    std::string atom,potential,mesh;
    int angular,Ne,order;
    double r0, rN,lambda = 0.0;
    int flag = read_data(argv,potential,atom,angular,Ne, r0,rN,order,mesh,lambda);
    std::cout<<"Your initial data:"<<std::endl;
    std::cout<<"Potential: "<<potential<<std::endl;
    if(potential=="Atomic-Screened"){
        std::cout<<"Lambda value: "<<lambda<<std::endl;
    }
    std::cout<<"Atom: "<<atom<<std::endl;
    std::cout<<"Angular Momentum (l): "<<angular<<std::endl;
    std::cout<<"Mesh: "<<mesh<<std::endl;






}