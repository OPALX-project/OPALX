/*!
 * @file error.cpp
 * @author Matthias Frey
 * @date 25. - 26. October 2016
 * @details Compares single-level solve with multi-level solve using
 * the PoissonProblems class. Each run appends the \f$ L_{2}\f$-error
 * to the file \e l2_error.dat.
 * @brief Compares the single-level solve with the multi-level solve
 */

#include <fstream>
#include <iostream>
#include <sstream>

#include "PoissonProblems.h"

int main(int argc, char* argv[]) {
    
    std::stringstream cmd;
    cmd << "mpirun -np [#cores] testError [gridx] [gridy] [gridz] "
        << "[max. grid] [#levels] [hor. domain, e.g. 0 1] [vert. domain] "
        << "[long. domain] [NOPARTICLES/UNIFORM/GAUSSIAN/REAL]";
    
    if ( argc < 13 ) {
        std::cerr << cmd.str() << std::endl;
        return -1;
    }
    
    Ippl ippl(argc, argv);
    
    BoxLib::Initialize(argc, argv, false);
    
    int nr[BL_SPACEDIM] = {
        std::atoi(argv[1]),
        std::atoi(argv[2]),
        std::atoi(argv[3])
    };
    
    
    int maxGridSize = std::atoi(argv[4]);
    int nLevels = std::atoi(argv[5]);
    
    double l2error = 0.0;
    bool solved = true;
    std::vector<double> lower = { std::atof(argv[6]),
                                  std::atof(argv[8]),
                                  std::atof(argv[10])
    };
    
    std::vector<double> upper = { std::atof(argv[7]),
                                  std::atof(argv[9]),
                                  std::atof(argv[11])
    };
    
    PoissonProblems pp(nr, maxGridSize, nLevels, lower, upper);
    
    if ( std::strcmp(argv[12], "NOPARTICLES") == 0 )
        l2error = pp.doSolveNoParticles();
    else if ( std::strcmp(argv[12], "UNIFORM") == 0 )
        l2error = pp.doSolveParticlesUniform();
    else if ( std::strcmp(argv[12], "GAUSSIAN") == 0 ) {
        
        if ( argc != 16 ) {
            std::cerr << cmd.str() << " [#particles] [mean] [stddev]" << std::endl;
            return -1;
        }
        
        int nParticles = std::atoi(argv[13]);
        double mean = std::atof(argv[14]);
        double stddev = std::atof(argv[15]);    // standard deviation
        
        l2error = pp.doSolveParticlesGaussian(nParticles, mean, stddev);
    } else if ( std::strcmp(argv[12], "REAL") == 0 ) {
        
        if ( argc != 15 ) {
            std::cerr << cmd.str() << " [step] [REAL: H5 file]" << std::endl;
            return -1;
        }
        
        int step = std::atoi(argv[13]);
        
        l2error = pp.doSolveParticlesReal(step, argv[14]);
    } else {
        if ( ParallelDescriptor::MyProc() == 0 )
            std::cerr << "Not supported solver." << std::endl
                      << "Use:" << std::endl
                      << "- NOPARTICLES: rhs = -1 everywhere" << std::endl
                      << "- UNIFORM: As NOPARTICLES but with particles (single-core)"
                      << std::endl
                      << "- GAUSSIAN: rhs using particles" << std::endl;
        solved = false;
    }
    
    
    
    ParallelDescriptor::Barrier();
    
    if ( ParallelDescriptor::MyProc() == 0 && solved) {
        std::ofstream out("l2_error_" + std::string(argv[6]) + ".dat", std::ios::app);
        out << nLevels << " " << l2error << std::endl;
        out.close();
    }
    
//     BoxLib::Finalize();
    
    return 0;
}