/*!
 * @file toFile.cpp
 * @author Matthias Frey
 * @date 3. Jan. 2017
 * @brief Read in a particle distribution from a
 *        H5 file and dump it ot a text file.
 * @details It dumps a file in the format of chapter
 *          11.3 in the OPAL manual.
 */
#include "Distribution.h"
#include "Ippl.h"

int main(int argc, char* argv[]) {
    
    Ippl ippl(argc, argv);
    
    std::string file = argv[1];
    int step = std::atoi(argv[2]);
    
    
    Distribution distr;
    
    distr.readH5(file, step);
    
    distr.print2file("distribution.dat");
    
    return 0;
}