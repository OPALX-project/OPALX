/* Author: Matthias Frey
 * 25. - 26. October 2016
 * 
 * Solve \Lap\ph = \rho
 * on the domain [0, 1] x [0, 1] x [0, 1]
 * for different ways of \rho initializations
 * 
 * doSolveNoParticles:          \rho = -1 everywhere (no particles)
 * doSolveParticlesUniform:     \rho = -1 everywhere (initialized by particles, idea: we generate particles on the
 *                              cell center and then do an AssignDensity-function call --> \rho = -1 everywhere
 * doSolveParticlesGaussian:    \rho is a Gaussian distribution initialized by particles
 * doSolveParticlesReal:        Read in a H5 cyclotron file and use its particle distribution for \rho
 */

#include <fstream>
#include <iostream>

#include "PoissonProblems.h"

int main(int argc, char* argv[]) {
    
    if ( argc != 8 ) {
        std::cerr << "mpirun -np [#cores] testError [gridx] [gridy] [gridz] "
                  << "[max. grid] [#levels] "
                  << "[NOPARTICLES/UNIFORM/GAUSSIAN/REAL] "
                  << "[#particles] [REAL: H5 file]" << std::endl;
        return -1;
    }
    
    BoxLib::Initialize(argc, argv, false);
    
    int nr[BL_SPACEDIM] = {
        std::atoi(argv[1]),
        std::atoi(argv[2]),
        std::atoi(argv[3])
    };
    
    
    int maxGridSize = std::atoi(argv[4]);
    int nLevels = std::atoi(argv[5]);
    int nParticles = std::atoi(argv[7]);
    
    double l2error = 0.0;
    bool solved = true;
    
    PoissonProblems pp(nr, maxGridSize, nLevels);
    
    if ( std::strcmp(argv[6], "NOPARTICLES") == 0 )
        l2error = pp.doSolveNoParticles();
    else if ( std::strcmp(argv[6], "UNIFORM") == 0 )
        l2error = pp.doSolveParticlesUniform();
    else if ( std::strcmp(argv[6], "GAUSSIAN") == 0 )
        l2error = pp.doSolveParticlesGaussian(nParticles);
    else if ( std::strcmp(argv[6], "REAL") == 0 )
        l2error = pp.doSolveParticlesReal(nParticles /* step */, argv[8]);
    else {
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
        std::ofstream out("l2_error.dat", std::ios::app);
        out << nLevels << " " << l2error << std::endl;
        out.close();
    }
    
    BoxLib::Finalize();
    
    return 0;
}