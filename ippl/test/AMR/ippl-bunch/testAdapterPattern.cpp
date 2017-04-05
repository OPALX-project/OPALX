#include "Ippl.h"

#include "PartBunch.h"
#include "AmrPartBunch.h"
#include "BoxLibClasses.h"

#include <memory>


void createRandomParticles(PartBunchBase<double, 3>* pbase, int N, int myNode, int seed = 1) {
    srand(seed);
    for (int i = 0; i < N; ++i) {
        pbase->createWithID(myNode * N + i + 1);
        pbase->Q[i] = (double)rand() / RAND_MAX;
    
        pbase->R[i][0] = (double)rand() / RAND_MAX;
        pbase->R[i][1] = (double)rand() / RAND_MAX;
        pbase->R[i][2] = (double)rand() / RAND_MAX;  
    }
}

int main(int argc, char** argv) {
    
    Ippl ippl(argc, argv);
    
    std::cout << "Init PartBunch" << std::endl
              << "--------------" << std::endl;
    
    std::unique_ptr<PartBunchBase<double, 3> > bunch(
        new PartBunch(
            new IpplParticleBase<ParticleSpatialLayout<double, 3> >(
                new ParticleSpatialLayout<double, 3>()
            )
        )
    );
    
    
    createRandomParticles(bunch.get(), 10, Ippl::myNode());
    
    std::cout << bunch->getTotalNum() << " " << bunch->getLocalNum() << std::endl;
    
    bunch->update();
    
    std::cout << bunch->getTotalNum() << " " << bunch->getLocalNum() << std::endl;
    
    std::cout << "bunch->R[0] = " << bunch->R[0] << std::endl;
    
    double two = 2.0;
    bunch->R *= Vector_t(two);
    
    std::cout << "2 * bunch->R[0] = " << bunch->R[0] << std::endl;
    std::cout << "max(bunch->R) = " << max(bunch->R) << std::endl;
    
    std::cout << std::endl << "Init AmrPartBunch" << std::endl
              << "-----------------" << std::endl;
    
    bunch.reset(
        new AmrPartBunch(
            new BoxLibParticle<BoxLibLayout<double, 3> >(
                new BoxLibLayout<double, 3>()
            )
        )
    );
    
    createRandomParticles(bunch.get(), 10, Ippl::myNode());
    
    
    std::cout << bunch->getTotalNum() << " " << bunch->getLocalNum() << std::endl;
    
    bunch->update();
    
    std::cout << bunch->getTotalNum() << " " << bunch->getLocalNum() << std::endl;
    
    std::cout << "bunch->R[0] = " << bunch->R[0] << std::endl;
    
    bunch->R *= Vector_t(two);
    
    std::cout << "2 * bunch->R[0] = " << bunch->R[0] << std::endl;
    std::cout << "max(bunch->R) = " << max(bunch->R) << std::endl;
    
    
    return 0;
}