#include "UserInputData.hpp"

//Reads the user input Data:
//FEM Methd: Fixed Elements
//Atomic Model: Confined
//Kind Of Confinement: Hard-Walls
//Rc: 1.0
//Atom: He;
//Charge +1;
//Angular Momentum: 0
//*** ONLY FOR DEVELOPERS****
//RInfty: 200.0;
//Integration: Numeric;
/*int read_data(char *input_name[], std::string &femModel, 
int Ne, int order, std::string &atomicModel, 
std::string &confType, double &rC, double &wallValue, double &lambda, int &charge, std::string &atomSym, int &angular, double &rInfty, std::string &mesh, std::string &integrals)*/
int read_data(char *input_name[], std::string &femModel, 
int &Ne, int &order, std::string &atomicModel, 
std::string &confType, double &rC, 
double &wallValue, double &lambda, 
int &charge, std::string &atomSym,std::string &multi, int &angular, 
double r0, double &rInfty, 
std::string &mesh, std::string &integrals)
{
    angular = -1; Ne = 0; order = 0; r0 = 0.0;
    std::ifstream inputfile;
    std::string filename="input.dat";
    inputfile.open(input_name[1]);
    

    if(!inputfile.is_open())
    {
        return 1;
    }
    std::string line;
    getline(inputfile, line,':');
    inputfile>>femModel;

    if(femModel=="Fixed-Elements"){
        inputfile>>std::ws;
        getline(inputfile, line,':');
        inputfile>>Ne;
    }
    inputfile>>std::ws;
    getline(inputfile, line,':');
    inputfile>>order;

    inputfile>>std::ws;
    getline(inputfile,line,':');
    inputfile >> atomicModel;
    if(atomicModel=="Confined")
    {
        inputfile>>std::ws;
        getline(inputfile,line,':');
        inputfile >> confType;
        if(confType=="Soft-Walls"){
            inputfile>>std::ws;
            getline(inputfile,line,':');
            inputfile>>wallValue;

            inputfile>>std::ws;
            getline(inputfile,line,':');
            inputfile>>rC;

        }
        else{
            inputfile>>std::ws;
            getline(inputfile,line,':');
            inputfile>>rC;

        }
    }
    else if(atomicModel=="Plasma"){
        inputfile>>std::ws;
        getline(inputfile,line,':');
        inputfile>>lambda;
    }
    inputfile>>std::ws;
    getline(inputfile,line,':');
    inputfile >> atomSym;

    inputfile>>std::ws;
    getline(inputfile,line,':');
    inputfile >> multi;


    inputfile>>std::ws;
    getline(inputfile,line,':');
    inputfile >> charge;

    inputfile>>std::ws;
    getline(inputfile,line,':');
    inputfile >> angular;

    inputfile>>std::ws;
    getline(inputfile,line,':');
    inputfile >> rInfty;

    inputfile>>std::ws;
    getline(inputfile,line,':');
    inputfile >> mesh;

    if(femModel=="Fixed Points"){

        inputfile>>std::ws;
        getline(inputfile,line,':');
        inputfile >> integrals;

    }


    inputfile.close();

    return 0;
}  

int validateUsrData(char *input_name[], std::string &femModel, 
    int &Ne, int &order, std::string &atomicModel, 
    std::string &confType, double &rC, 
    double &wallValue, double &lambda, 
    int &charge, std::string &atomSym, int &angular, 
    double r0, double &rInfty, 
    std::string &mesh, std::string &integrals){
    if(femModel!="Fixed-Elements" || femModel!="Fixed-Points"){
        std::cout<<"ERROR: FEM Model not recognized\n";
        return 1;
    }
    else if(Ne==0 || order ==0){
        std::cout<<"ERROR: Number of elements or interpolation order is Zero\n";
        return 1;
    }
    


    return 0;

}
void printUsrData(std::string femModel, 
    int Ne, int order, std::string atomicModel, 
    std::string confType, double rC, 
    double wallValue, double lambda, 
    int charge, std::string atomSym,std::string multi, int angular, 
    double r0, double rN, 
    std::string mesh, std::string integrals){
    printf("**************************************************\n");
    std::cout<<"Your initial data:"<<std::endl;
    std::cout<<"FEM Model: "<<femModel<<std::endl;
    std::cout<<"Number of Elements: "<<Ne<<std::endl;
    std::cout<<"Number of Points: "<<Ne*order + 1<<std::endl;
    std::cout<<"Order: "<<order<<std::endl;
    std::cout<<"Atomic Model: "<<atomicModel<<std::endl;
    if(atomicModel=="Confined"){
        if(confType=="Soft-Walls"){
            std::cout<<"Confinement: "<<confType<<std::endl;
            std::cout<<"Wall Value: "<<wallValue<<std::endl;
            std::cout<<"Rc: "<<rC<<std::endl;
        }
        else{
            std::cout<<"Confinement: "<<confType<<std::endl;
            std::cout<<"Rc: "<<rC<<std::endl;
        }
            
    }
    else if(atomicModel=="Plasma"){
        std::cout<<"Lambda Value: "<<lambda<<std::endl;
            
    }
    std::cout<<"Atom: "<<atomSym<<std::endl;
    std::cout<<"Charge: "<<charge<<std::endl;
    std::cout<<"Multiplicity: "<<multi<<std::endl;
    std::cout<<"Angular Momentum: "<<angular<<std::endl;
    std::cout<<"rInfty: "<<rN<<std::endl;

    std::cout<<"Mesh: "<<mesh<<std::endl;
    printf("**************************************************\n");

 
}