#include <iostream>
#include <fstream>
#if !defined(_USER_INPUT_H_)
#define _USER_INPUT_H_
#include <iostream>
#include <fstream>
int read_data(char *input_name[],std::string &potential,std::string &atom, int &angular,int &Ne, double &r0, double &rc,int &order, std::string &mesh,double &lambda);

#endif // _USER_INPUT_H_
