#include <iostream>
#include <fstream>
#if !defined(_USER_INPUT_H_)
#define _USER_INPUT_H_
#include <iostream>
#include <fstream>
int read_data(char *input_name[], std::string &femModel, 
int &Ne, int &order, std::string &atomicModel, 
std::string &confType, double &rC, 
double &wallValue, double &lambda, 
int &charge, std::string &atomSym,std::string &multi, int &angular, 
double r0, double &rInfty, 
std::string &mesh, std::string &integrals);

int validateUsrData(char *input_name[], std::string &femModel, 
int &Ne, int &order, std::string &atomicModel, 
std::string &confType, double &rC, 
double &wallValue, double &lambda, 
int &charge, std::string &atomSym, int &angular, 
double r0, double &rInfty, 
std::string &mesh, std::string &integrals);

void printUsrData(std::string femModel, 
int Ne, int order, std::string atomicModel, 
std::string confType, double rC, 
double wallValue, double lambda, 
int charge, std::string atomSym,std::string multi, int angular, 
double r0, double rN, 
std::string mesh, std::string integrals);
#endif // _USER_INPUT_H_
