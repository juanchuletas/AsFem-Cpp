#include "searchelements.hpp"


int getAtomicNumber(std::string target){

    int totElements = (int) (sizeof(chem_elements)/sizeof(chem_elements[0]));
    for(int z=1; z<totElements; z++){

        if(chem_elements[z]==target){
            return z;
        }

    }

    std::cout<<"Chemical Element not found\n";
    return 0;

}