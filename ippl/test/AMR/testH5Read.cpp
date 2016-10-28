#include <iostream>

#include "Distribution.h"

int main(int argc, char* argv[]) {
    
    Ippl ippl(argc, argv);
    
    Distribution distr;
    
    
    std::string file = "../../../../../finalsim.h5";
    distr.readH5(file, 0);
    
    return 0;
}